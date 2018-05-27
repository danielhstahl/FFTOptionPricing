#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include "RungeKutta.h"
#include "OptionPricing.h"
#include <chrono>
#include "fft.h"
#include <deque>
#include <iostream>
#include "CharacteristicFunctions.h"
#include "cuckoo.h"
#include "utils.h"
#include "OptionCalibration.h"


const auto numSims=1500;
const auto nestSize=25;
const double tol=.000001;
template<typename CF>
void printResults(
    const std::string& nameToWrite,  
    CF&& cf, 
    const std::vector<double>& strikes,
    std::vector<double>& objParms,
    const std::vector<swarm_utils::upper_lower<double>>& constraints,
    const std::vector<std::string>& paramNames,
    double r,
    double T, 
    double stock, 
    double strikeMultiplier,
    std::vector<double>&& uArray
){
    if(constraints.size()!=paramNames.size()){
        std::cout<<"REQUIRES constraints and paramNames SAME SIZE"<<std::endl;
    }
    if(objParms.size()!=paramNames.size()){
        std::cout<<"REQUIRES objParms and paramNames SAME SIZE"<<std::endl;
    }
    double discount=exp(-r*T);
    double maxStrike=strikes.back()*strikeMultiplier;
    double minStrike=stock/maxStrike;
    int n=strikes.size();
    std::vector<double> KArray(n+2);
    for(int i=0; i<n; ++i){
        KArray[i+1]=strikes[n-i-1];
    }
    KArray[0]=stock*50; //endpoints for creating option price
    KArray[n+1]=stock*.03;
    int numU=256;
    auto optionPrices=optionprice::FangOostCallPrice(stock, KArray, r, T, numU, [&](const auto& u){
        return exp(r*T*u+cf(objParms)(u));
    });
    auto getReducedAndReversed=[](const std::vector<double>& arr){  
        auto arrTmp=std::vector<double>(arr.begin()+1, arr.end()-1);
        std::reverse(std::begin(arrTmp), std::end(arrTmp));
        return arrTmp;
    };
    auto observedO=getReducedAndReversed(optionPrices);
    if(observedO.size()!=strikes.size()){
        std::cout<<"SIZE OF OPTIONS (observedO, "<<observedO.size()<<") MUST EQUAL STRIKES (strikes, "<<strikes.size()<<")"<<std::endl;
    }
    auto s=optioncal::getOptionSpline(strikes, observedO, stock, discount, minStrike, maxStrike);
    auto valOrZero=[](const auto& v){
        return v>0.0?v:0.0;
    };
    
    double minLogStrike=log(minStrike);
    double maxLogStrike=log(maxStrike/stock);
    int m=256;
    double logdk=(maxLogStrike-minLogStrike)/(double)(m-1);
    std::vector<double> dkLogArray(m+2);
    dkLogArray=futilities::for_each_parallel_subset(std::move(dkLogArray), 1, 1, [&](const auto& v, const auto& index){
        return exp(maxLogStrike-(index-1)*logdk);
    });
    dkLogArray[0]=50;
    dkLogArray[m+1]=.03;

    auto optionPricesLogDK=optionprice::FangOostCallPrice(1.0, dkLogArray, r, T, numU, [&](const auto& u){
        return exp(r*T*u+cf(objParms)(u));
    });
    std::ofstream splineFileCall;

    splineFileCall.open("./docs/calibration/"+nameToWrite+"Spline.csv");

    splineFileCall<<"strike, price, actual"<<std::endl;
    for(int i=0; i<m;++i){
        double logk=minLogStrike+i*logdk;
        double offset=1-exp(logk)*discount;
        splineFileCall<<logk-r*T<<", "<<valOrZero(s(exp(logk)))<<", "<<optionPricesLogDK[m-i]-valOrZero(offset)<<std::endl;
    }
    splineFileCall.close();
    
    int N=4096;
    const auto estimateOfPhi=optioncal::generateFOEstimate(strikes, observedO, stock, r, T, minStrike, maxStrike);

    auto phis=estimateOfPhi(N, uArray);
    std::ofstream numericIntegration;
    std::ofstream optimalParameters;
    numericIntegration.open("./docs/calibration/"+nameToWrite+"Integration.csv");
    numericIntegration<<"u, estimateReal, estimateIm, exactReal, exactIm"<<std::endl;
    auto cfh=cf(objParms);
    for(int i=0; i<phis.size(); ++i){ 
        const auto u=uArray[i];
        const auto actualPhi=cfh(std::complex<double>(1.0, u));
        numericIntegration<<u<<", "<<phis[i].real()<<", "<<phis[i].imag()<<", "<<actualPhi.real()<<","<<actualPhi.imag()<<std::endl;
    }
    numericIntegration.close();
   
    auto objFn=optioncal::getObjFn_arr(
        std::move(phis),
        std::move(cf),
        std::move(uArray)
    );
    
    const auto results=cuckoo::optimize(objFn, constraints, nestSize, numSims, tol, 42);

    auto params=std::get<swarm_utils::optparms>(results);
    auto objFnRes=std::get<swarm_utils::fnval>(results);

    optimalParameters.open("./docs/calibration/"+nameToWrite+"Parameters.csv");
    optimalParameters<<"paramater, actual, optimal"<<std::endl;
    int numParameters=constraints.size();
    for(int i=0; i<numParameters;++i){
        optimalParameters<<paramNames[i]<<", "<<objParms[i]<<", "<<params[i]<<std::endl;
    }
    optimalParameters.close();
}

template<typename CF>
void printResults(
    const std::string& nameToWrite,  
    CF&& cf, 
    const std::vector<double>& strikes,
    std::vector<double>& objParms,
    const std::vector<swarm_utils::upper_lower<double>>& constraints,
    const std::vector<std::string>& paramNames,
    double r,
    double T, 
    double stock,
    double strikeMultiplier
){
    //works well on 0 to 2pi
    auto getU=[](const auto& N){
        double du= 2.0*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    auto uArray=getU(15);
    printResults(nameToWrite, std::move(cf), strikes, objParms, constraints,
    paramNames, r, T, stock, strikeMultiplier, std::move(uArray));
}

TEST_CASE("Test BS Large U", "[OptionCalibration]"){
    double stock=10.0;
    double r=0;
    double sigma=.3;
    double T=1.0;
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
    };
    std::vector<std::string> paramNames={
        "sigma"
    };
    std::vector<double> strikes={4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    auto cf=[T](const auto& objParams){
        return [=](const auto&u){
            auto sigma=objParams[0];
            return (-sigma*sigma*u*.5+sigma*sigma*u*u*.5)*T;
        };
    };
    std::vector<double> objParams={sigma};
    std::vector<double> uArray={-20, -15, -10, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 10, 15, 20};
    printResults("BlackScholesU", std::move(cf), strikes, objParams, constraints, paramNames, r, T, stock, 3.0, std::move(uArray));
}


TEST_CASE("Test BS", "[OptionCalibration]"){
    double stock=10.0;
    double r=.05;
    double sigma=.3;
    double T=1.0;
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
    };
    std::vector<std::string> paramNames={
        "sigma"
    };
    std::vector<double> strikes={4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    auto cf=[T](const auto& objParams){
        auto sigma=objParams[0];
        return [=](const auto&u){
            return (-sigma*sigma*u*.5+sigma*sigma*u*u*.5)*T;
        };
    };
    std::vector<double> objParams={sigma};
    printResults("BlackScholes", std::move(cf), strikes, objParams, constraints, paramNames, r, T, stock, 3.0);
}

TEST_CASE("Test Heston", "[OptionCalibration]"){
    double r=0;
    double T=1.0;
    double b=.0398;
    double a=1.5768;
    double c=.5751;
    double rho=-.5711;
    double v0=.0175;
    //convert to extended CGMY
    auto sig=sqrt(b);
    auto speed=a;
    auto stock=178.46;  
    //auto kappa=speed;//long run tau of 1
    auto v0Hat=v0/b;
    auto adaV=c/sqrt(b);
    std::vector<double> objParms={sig, speed, adaV, rho, v0Hat};

    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 4.0),
        swarm_utils::upper_lower<double>(-1.0, 1.0),
        swarm_utils::upper_lower<double>(0, 2),
    };
    std::vector<std::string> paramNames={
        "sigma","speed", "adaV", "rho", "v0Hat"
    };

    std::vector<double> strikes={
        95,100,130,150,160,165,170,175,185,190,195,200,210,240,250
    };    
    auto cf=[T](
        const std::vector<double>& params
    ){
        auto sig=params[0];
        auto speed=params[1];
        auto adaV=params[2];
        auto rho=params[3];
        auto v0Hat=params[4];
        return [=](const auto& u){
            return chfunctions::cirLogMGF(
                -chfunctions::cgmyLogRNCF(u, 0.0, 1.0, 1.0, .5, 0.0, sig),
                speed, 
                speed-adaV*rho*u*sig,
                adaV,
                T, 
                v0Hat
            );
        };
        
    };
    printResults("Heston", std::move(cf), strikes, objParms, constraints, paramNames, r, T, stock, 3.0);
}

TEST_CASE("Test CGMY", "[OptionCalibration]"){
    auto r=.04;
    auto sig=0.1;
    auto T=1.0;
    auto S0=178.0;  
    auto C=.1;
    auto G=2.4;
    auto M=10.7;
    auto Y=.5;

    std::vector<double> objParms={sig, C, G, M, Y};

    auto CFBase=[T](const std::vector<double>& params){
        auto sig=params[0];
        auto C=params[1];
        auto G=params[2];
        auto M=params[3];
        auto Y=params[4];
        return [=](const auto& u){
            return chfunctions::cgmyLogRNCF(u, C, G, M, Y, 0.0, sig)*T;
        };
        
    };
    std::vector<double> strikes={
        95,100,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<std::string> paramNames={
        "sig", "C", "G", "M", "Y"
    };
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 15),
        swarm_utils::upper_lower<double>(0.0, 15),
        swarm_utils::upper_lower<double>(-3, 2)
    };
    printResults("CGMY", std::move(CFBase), strikes, objParms, constraints, paramNames, r, T, S0, 10.0);
   
}

TEST_CASE("Test Merton", "[OptionCalibration]"){
    auto r=.04;//seems high
    auto sig=0.2;
    auto T=.25;
    auto S0=178.0;  
    auto lambda=.5;
    auto muJ=.05;
    auto sigJ=.05;
    double discount=exp(-r*T);

    std::vector<double> objParms={sig, lambda, muJ, sigJ};

    auto CFBase=[T](
        const std::vector<double>& params
    ){
        auto sig=params[0];
        auto lambda=params[1];
        auto muJ=params[2];
        auto sigJ=params[3];
        return [=](const auto& u){
            return chfunctions::mertonLogRNCF(u, lambda, muJ, sigJ, 0.0, sig)*T;
        };
        
    };
    std::vector<double> strikes={
        95,100,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<std::string> paramNames={
        "sig", "lambda", "muJ", "sigJ"
    };
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 1.5),
        swarm_utils::upper_lower<double>(0.0, 1.5)
    };
    printResults("Merton", std::move(CFBase), strikes, objParms, constraints, paramNames, r, T, S0, 10.0);
   
}
TEST_CASE("Test Merton with stochastic vol", "[OptionCalibration]"){
    auto r=.04;//seems high
    auto sig=0.2;
    auto T=.25;
    auto S0=178.0;  
    auto lambda=.5;
    auto muJ=.05;
    auto sigJ=.05;
    auto rho=-.5;
    auto speed=.5;
    auto adaV=.2;
    auto v0=.9;
    double discount=exp(-r*T);

    std::vector<double> objParms={sig, lambda, muJ, sigJ, speed, adaV, rho, v0};

    auto CFBase=[T](
        const std::vector<double>& params
    ){
        auto sig=params[0];
        auto lambda=params[1];
        auto muJ=params[2];
        auto sigJ=params[3];
        auto speed=params[4];
        auto adaV=params[5];
        auto rho=params[6];
        auto v0Hat=params[7];
        return [=](const auto& u){
            return chfunctions::cirLogMGF(
                -chfunctions::mertonLogRNCF(u, lambda, muJ, sigJ, 0.0, sig),
                speed, 
                speed-adaV*rho*u*sig,
                adaV,
                T, 
                v0Hat
            );
        };
        
    };
    std::vector<double> strikes={
        95,100,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<std::string> paramNames={
        "sig", "lambda", "muJ", "sigJ", "speed", "adaV", "rho", "v0"
    };
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 1.5),
        swarm_utils::upper_lower<double>(0.0, 1.5),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 4.0),
        swarm_utils::upper_lower<double>(-1.0, 1.0),
        swarm_utils::upper_lower<double>(0, 2)
    };
    printResults("MertonCorr", std::move(CFBase), strikes, objParms, constraints, paramNames, r, T, S0, 10.0);
   
}
auto cfLogGeneric(
    double T
){
    return [=](
        double lambda, 
        double muJ, double sigJ,
        double sigma, double v0, 
        double speed, double adaV, 
        double rho, double delta
    ){

        auto numODE=64;//hopefully this is sufficient
        //double speedTmp=speed;//copy here in order to move to move
        //const T& rho, const T& K, const T& H, const T& l
        auto alpha=chfunctions::AlphaOrBeta_move(0.0, speed, 0.0, 0.0);
        return [=, alpha=std::move(alpha), numODE=std::move(numODE)](const auto& u){
            //const T& rho, const T& K, const T& H, const T& l
            auto beta=chfunctions::AlphaOrBeta_move(
                -chfunctions::mertonLogRNCF(u, lambda, muJ, sigJ, 0.0, sigma), 
                -(speed+delta*lambda-u*rho*sigma*adaV),
                adaV*adaV, 
                lambda
            );
            
            auto expCF=chfunctions::exponentialCFBeta(u, delta);
            return chfunctions::logAffine(
                rungekutta::compute_efficient_2d(
                    T, numODE, 
                    std::vector<std::complex<double> >({0, 0}),
                    [&](
                        double t, 
                        const std::complex<double>& x1,  
                        const std::complex<double>& x2
                    ){
                        auto cfPart=(chfunctions::exponentialCFBeta(
                            x1, 
                            delta
                        )-1.0)*chfunctions::gaussCF(u, muJ, sigJ);
                        return std::vector<std::complex<double> >({
                            beta(x1, cfPart),
                            alpha(x1, cfPart)
                        });
                    }
                ),
                v0
            );
        };
    };
}
TEST_CASE("Test time changed merton with volatility jumps", "[OptionCalibration]"){
    auto r=.04;
    auto sig=0.1;
    auto T=1;
    auto S0=178;  
    auto lambda=.3;
    auto muJ=-.05;
    auto sigJ=.1;
    auto speed=.3;
    auto v0=.9;
    auto adaV=.2;
    auto rho=-.5;
    auto delta=.2;
    double discount=exp(-r*T);

    std::vector<double> objParms={lambda, muJ, sigJ, sig, v0, speed, adaV, rho, delta};

    auto CFBase=[T](
        const std::vector<double>& params
    ){
        auto lambda=params[0];
        auto muJ=params[1];
        auto sigJ=params[2];
        auto sigma=params[3];
        auto v0=params[4];
        auto speed=params[5];
        auto adaV=params[6];
        auto rho=params[7];
        auto delta=params[8];
        return cfLogGeneric(T)(lambda, muJ, sigJ, sigma, v0, speed, adaV, rho, delta);
    };
    std::vector<double> strikes={
        95,100,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<std::string> paramNames={
       "lambda", "muJ", "sigJ", "sigma", "v0", "speed", "adaV", "rho", "delta"
    };
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(-1.0, 1.0),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 1.0),
        swarm_utils::upper_lower<double>(0.2, 1.8),
        swarm_utils::upper_lower<double>(0.0, 3.0),
        swarm_utils::upper_lower<double>(0.0, 3.0),
        swarm_utils::upper_lower<double>(-1.0, 1.0),
        swarm_utils::upper_lower<double>(0, 2)
    };
    printResults("MertonJumps", std::move(CFBase), strikes, objParms, constraints, paramNames, r, T, S0, 20.0);
   
}

TEST_CASE("Test spline for actual data", "[OptionCalibration]"){
    std::vector<double> strikes={
        95,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<double> optionPrices={
        85,51.5,35.38,28.3,25.2,22.27,19.45,14.77,12.75,11,9.35,6.9,2.55,1.88
    }; 
    auto stock=178.46;
    auto r=.003;
    auto T=1.0;
    double discount=exp(-r*T);
    double maxStrike=strikes.back()*4.5;
    double minStrike=stock/maxStrike;
    auto s=optioncal::getOptionSpline(strikes, optionPrices, stock, discount, minStrike, maxStrike);
    std::ofstream splineRealData;
    std::ofstream splineActualData;
    auto valOrZero=[](const auto& v){
        return v>0.0?v:0.0;
    };
    splineActualData.open("./docs/calibration/ActualDataSpline.csv");
    splineActualData<<"strike, price"<<std::endl;
    for(int i=0; i<strikes.size();++i){
        double offset=1-strikes[i]*discount/stock;
        splineActualData<<log(strikes[i]/stock)-r*T<<","<<optionPrices[i]/stock-valOrZero(offset)<<std::endl;
    }
    splineActualData.close();

    splineRealData.open("./docs/calibration/RealDataSpline.csv");
    splineRealData<<"strike, price"<<std::endl;
    int m=256;
    
    double minLogStrike=log(minStrike);
    double maxLogStrike=log(maxStrike/stock);
    
    double logdk=(maxLogStrike-minLogStrike)/(double)(m-1);
    for(int i=0; i<m;++i){
        double logk=minLogStrike+i*logdk;
        double offset=1-exp(logk)*discount;
        splineRealData<<logk-r*T<<", "<<valOrZero(s(exp(logk)))<<std::endl;
    }
    splineRealData.close();
}

TEST_CASE("Heston on actual data", "[OptionCalibration]"){
    std::vector<double> strikes={
        95,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<double> optionPrices={
        85,51.5,35.38,28.3,25.2,22.27,19.45,14.77,12.75,11,9.35,6.9,2.55,1.88
    }; 
    auto stock=178.46;
    auto r=.003;
    auto T=1.0;
    double maxStrike=strikes.back()*10.0;
    double minStrike=stock/maxStrike;
    auto foEstimate=optioncal::generateFOEstimate(strikes, optionPrices, stock, r, T, minStrike, maxStrike);
    
    std::ofstream optimalParameters;
    std::ofstream optimalValue;
    auto getU=[](const auto& N){
        double du= 2.0*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    auto uArray=getU(15);
    int N=1024;
    auto phis=foEstimate(N, uArray);
   
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 4.0),
        swarm_utils::upper_lower<double>(-1.0, 1.0),
        swarm_utils::upper_lower<double>(0, 2),
    };
    std::vector<std::string> paramNames={
        "sigma","speed", "adaV", "rho", "v0Hat"
    };
  
    auto cf=[T](
        
        const std::vector<double>& params
    ){
        auto sig=params[0];
        auto speed=params[1];
        auto adaV=params[2];
        auto rho=params[3];
        auto v0Hat=params[4];
        return [&](const auto& u){
            return chfunctions::cirLogMGF(
                -chfunctions::cgmyLogRNCF(u, 0.0, 1.0, 1.0, .5, 0.0, sig),
                speed, 
                speed-adaV*rho*u*sig,
                adaV,
                T, 
                v0Hat
            );
        };
        
    };
    auto objFn=optioncal::getObjFn_arr(
        std::move(phis),
        std::move(cf),
        std::move(uArray)
    );

    const auto results=cuckoo::optimize(objFn, constraints, nestSize, numSims, tol, 42);

    auto params=std::get<swarm_utils::optparms>(results);
    auto objFnRes=std::get<swarm_utils::fnval>(results);

    optimalParameters.open("./docs/calibration/HestonRealParameters.csv");
    optimalParameters<<"paramater, estimate "<<std::endl;
    int numParameters=constraints.size();
    for(int i=0; i<numParameters;++i){
        optimalParameters<<paramNames[i]<<", "<<params[i]<<std::endl;
    }
    optimalParameters.close();

    optimalValue.open("./docs/calibration/HestonMSE.csv");
    optimalValue<<"value"<<std::endl;
    optimalValue<<objFnRes<<std::endl;
    optimalValue.close();
}
TEST_CASE("CGMY on actual data", "[OptionCalibration]"){
    std::vector<double> strikes={
        95,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<double> optionPrices={
        85,51.5,35.38,28.3,25.2,22.27,19.45,14.77,12.75,11,9.35,6.9,2.55,1.88
    }; 
    auto stock=178.46;
    auto r=.003;
    auto T=1.0;
    double maxStrike=strikes.back()*10.0;
    double minStrike=stock/maxStrike;
    auto foEstimate=optioncal::generateFOEstimate(strikes, optionPrices, stock, r, T, minStrike, maxStrike);
    
    std::ofstream optimalParameters;
    std::ofstream optimalValue;
    auto getU=[](const auto& N){
        double du= 2.0*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    auto uArray=getU(15);
    int N=1024;
    auto phis=foEstimate(N, uArray);
   

    auto cf=[T](const std::vector<double>& params){
        auto sig=params[0];
        auto C=params[1];
        auto G=params[2];
        auto M=params[3];
        auto Y=params[4];
        return [=](const auto& u){
            return chfunctions::cgmyLogRNCF(u, C, G, M, Y, 0.0, sig)*T;
        };
    };
    
    std::vector<std::string> paramNames={
        "sig", "C", "G", "M", "Y"
    };
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 15),
        swarm_utils::upper_lower<double>(0.0, 15),
        swarm_utils::upper_lower<double>(-3, 2)
        
    };
    auto objFn=optioncal::getObjFn_arr(
        std::move(phis),
        std::move(cf),
        std::move(uArray)
    );

    const auto results=cuckoo::optimize(objFn, constraints, nestSize, numSims, tol, 42);

    auto params=std::get<swarm_utils::optparms>(results);
    auto objFnRes=std::get<swarm_utils::fnval>(results);

    optimalParameters.open("./docs/calibration/CGMYRealParameters.csv");
    optimalParameters<<"paramater, estimate "<<std::endl;
    int numParameters=constraints.size();
    for(int i=0; i<numParameters;++i){
        optimalParameters<<paramNames[i]<<", "<<params[i]<<std::endl;
    }
    optimalParameters.close();

    optimalValue.open("./docs/calibration/CGMYMSE.csv");
    optimalValue<<"value"<<std::endl;
    optimalValue<<objFnRes<<std::endl;
    optimalValue.close();
}

TEST_CASE("Time changed CGMY on actual data", "[OptionCalibration]"){
    std::vector<double> strikes={
        95,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<double> optionPrices={
        85,51.5,35.38,28.3,25.2,22.27,19.45,14.77,12.75,11,9.35,6.9,2.55,1.88
    }; 
    auto stock=178.46;
    auto r=.003;
    auto T=1.0;
    double maxStrike=strikes.back()*10.0;
    double minStrike=stock/maxStrike;
    auto foEstimate=optioncal::generateFOEstimate(strikes, optionPrices, stock, r, T, minStrike, maxStrike);
    
    std::ofstream optimalParameters;
    std::ofstream optimalValue;
    auto getU=[](const auto& N){
        double du= 2.0*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    auto uArray=getU(15);
    int N=1024;
    auto phis=foEstimate(N, uArray);
   

    auto cf=[T]( const std::vector<double>& params){
        auto sig=params[0];
        auto C=params[1];
        auto G=params[2];
        auto M=params[3];
        auto Y=params[4];
        auto speed=params[5];
        auto adaV=params[6];
        auto rho=params[7];
        auto v0=params[8];
        return [=](const auto& u){
            return chfunctions::cirLogMGF(
                -chfunctions::cgmyLogRNCF(u, C, G, M, Y, 0.0, sig),
                speed, 
                speed-adaV*rho*u*sig,
                adaV,
                T, 
                v0
            );
        };
        
    };
    
    std::vector<std::string> paramNames={
        "sig", "C", "G", "M", "Y", "speed", "adaV", "rho", "v0"
    };
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 15),
        swarm_utils::upper_lower<double>(0.0, 15),
        swarm_utils::upper_lower<double>(-3, 2),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 4.0),
        swarm_utils::upper_lower<double>(-1.0, 1.0),
        swarm_utils::upper_lower<double>(0, 2)
    };
    auto objFn=optioncal::getObjFn_arr(
        std::move(phis),
        std::move(cf),
        std::move(uArray)
    );

    const auto results=cuckoo::optimize(objFn, constraints, nestSize, numSims, tol, 42);

    auto params=std::get<swarm_utils::optparms>(results);
    auto objFnRes=std::get<swarm_utils::fnval>(results);

    optimalParameters.open("./docs/calibration/CGMYTimeChangeRealParameters.csv");
    optimalParameters<<"paramater, estimate "<<std::endl;
    int numParameters=constraints.size();
    for(int i=0; i<numParameters;++i){
        optimalParameters<<paramNames[i]<<", "<<params[i]<<std::endl;
    }
    optimalParameters.close();

    optimalValue.open("./docs/calibration/CGMYTimeChangeMSE.csv");
    optimalValue<<"value"<<std::endl;
    optimalValue<<objFnRes<<std::endl;
    optimalValue.close();
}

TEST_CASE("Time changed Merton", "[OptionCalibration]"){
    std::vector<double> strikes={
        95,130,150,160,165,170,175,185,190,195,200,210,240,250
    }; 
    std::vector<double> optionPrices={
        85,51.5,35.38,28.3,25.2,22.27,19.45,14.77,12.75,11,9.35,6.9,2.55,1.88
    }; 
    auto stock=178.46;
    auto r=.003;
    auto T=1.0;
    double maxStrike=strikes.back()*10.0;
    double minStrike=stock/maxStrike;
    auto foEstimate=optioncal::generateFOEstimate(strikes, optionPrices, stock, r, T, minStrike, maxStrike);
    
    std::ofstream optimalParameters;
    std::ofstream optimalValue;
    auto getU=[](const auto& N){
        double du= 2.0*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    auto uArray=getU(15);
    int N=1024;
    auto phis=foEstimate(N, uArray);
   

    auto cf=[T](
        
        const std::vector<double>& params
    ){
        auto sig=params[0];
        auto lambda=params[1];
        auto muJ=params[2];
        auto sigJ=params[3];
        auto speed=params[4];
        auto adaV=params[5];
        auto rho=params[6];
        auto v0=params[7];
        return [=](const auto& u){
            return chfunctions::cirLogMGF(
                -chfunctions::mertonLogRNCF(u, lambda, muJ, sigJ, 0.0, sig),
                speed, 
                speed-adaV*rho*u*sig,
                adaV,
                T, 
                v0
            );
        };
        
    };

    std::vector<std::string> paramNames={
        "sig", "lambda", "muJ", "sigJ", "speed", "adaV", "rho", "v0"
    };
    std::vector<swarm_utils::upper_lower<double>> constraints={
        swarm_utils::upper_lower<double>(0.0, .6),
        swarm_utils::upper_lower<double>(0.0, .5),
        swarm_utils::upper_lower<double>(0.0, 1.5),
        swarm_utils::upper_lower<double>(0.0, 1.5),
        swarm_utils::upper_lower<double>(0.0, 2.0),
        swarm_utils::upper_lower<double>(0.0, 4.0),
        swarm_utils::upper_lower<double>(-1.0, 1.0),
        swarm_utils::upper_lower<double>(0, 2)
    };
    
    auto objFn=optioncal::getObjFn_arr(
        std::move(phis),
        std::move(cf),
        std::move(uArray)
    );

    const auto results=cuckoo::optimize(objFn, constraints, nestSize, numSims, tol, 42);

    auto params=std::get<swarm_utils::optparms>(results);
    auto objFnRes=std::get<swarm_utils::fnval>(results);

    optimalParameters.open("./docs/calibration/MertonTimeChangeRealParameters.csv");
    optimalParameters<<"paramater, estimate "<<std::endl;
    int numParameters=constraints.size();
    for(int i=0; i<numParameters;++i){
        optimalParameters<<paramNames[i]<<", "<<params[i]<<std::endl;
    }
    optimalParameters.close();

    optimalValue.open("./docs/calibration/MertonTimeChangeMSE.csv");
    optimalValue<<"value"<<std::endl;
    optimalValue<<objFnRes<<std::endl;
    optimalValue.close();
}