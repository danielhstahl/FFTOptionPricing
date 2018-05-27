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

#include "OptionCalibration.h"
double BSCall(double S0, double discount, double k, double sigma, double T){
    double s=sqrt(2.0);
    double sigmasqrt=sqrt(T)*sigma;
    auto d1=log(S0/(discount*k))/(sigmasqrt)+sigmasqrt*.5;
    return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigmasqrt)/s)));
}
double BSPut(double S0, double discount, double k, double sigma, double T){
    return BSCall(S0, discount, k, sigma, T)-S0+k*discount;
}
double BSCallDelta(double S0, double discount, double k, double sigma, double T){
    double s=sqrt(2.0);
    auto d1=log(S0/(discount*k))/(sigma*sqrt(T))+sigma*.5*sqrt(T);
    return .5+.5*erf(d1/s);
}
double BSPutDelta(double S0, double discount, double k, double sigma, double T){
    return BSCallDelta(S0, discount, k, sigma, T)-1.0;
}
double BSGamma(double S0, double k, double sigma, double r, double T){
    double s=sqrt(2.0);
    auto d1=(log(S0/k)+r*T+sigma*sigma*.5*T)/(sigma*sqrt(T));
    return exp(-d1*d1/2.0)/(s*sqrt(M_PI)*sigma*sqrt(T)*S0);
}

double BSCallTheta(double S0, double k, double sigma, double r, double T){
    double s=sqrt(2.0);
    auto d2=(log(S0/k)+r*T-sigma*sigma*.5*T)/(sigma*sqrt(T));
    auto d1=d2+sigma*sqrt(T);
    auto dnorm=exp(-d1*d1/2.0)/(s*sqrt(M_PI));
    return -S0*dnorm*sigma/(2.0*sqrt(T))-r*k*exp(-r*T)*(.5+.5*erf(d2/s));
}
double BSPutTheta(double S0, double k, double sigma, double r, double T){
    return BSCallTheta(S0, k, sigma, r, T)+r*k*exp(-r*T);
}

TEST_CASE("FSTSCall", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    auto discount=exp(-r*T);
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto payoff=[&](const auto& assetValue){
        return assetValue>K?assetValue-K:0.0;
    };
    auto getAssetPrice=[&](const auto& logAsset){
        return K*exp(logAsset);
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FSTSPrice(numX, xmax, 
        r, T, 
        getAssetPrice,
        payoff,
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BSCall(myXDomain[i], discount, K, sig, T)).epsilon(.0001));
    }
}


TEST_CASE("FSTSCallDelta", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    auto discount=exp(-r*T);
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto payoff=[&](const auto& assetValue){
        return assetValue>K?assetValue-K:0.0;
    };
    auto getAssetPrice=[&](const auto& logAsset){
        return K*exp(logAsset);
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myDelta=optionprice::FSTSDelta(numX, xmax, 
        r, T, 
        getAssetPrice,
        payoff,
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.4;
    int mxX=(int)numX*.6;
    for(;i<mxX; ++i){
        REQUIRE(myDelta[i]==Approx(BSCallDelta(myXDomain[i], discount, K, sig, T)).epsilon(.0001));
    }
}

TEST_CASE("FSTSCallTheta", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    auto payoff=[&](const auto& assetValue){
        return assetValue>K?assetValue-K:0.0;
    };
    auto getAssetPrice=[&](const auto& logAsset){
        return K*exp(logAsset);
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myTheta=optionprice::FSTSTheta(numX, xmax, 
        r, T, 
        getAssetPrice,
        payoff,
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    
    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myTheta[i]==Approx(BSCallTheta(myXDomain[i], K, sig, r, T)).epsilon(.0001));
    }
}


TEST_CASE("FSTSCallGamma", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    auto payoff=[&](const auto& assetValue){
        return assetValue>K?assetValue-K:0.0;
    };
    auto getAssetPrice=[&](const auto& logAsset){
        return K*exp(logAsset);
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myGamma=optionprice::FSTSGamma(numX, xmax, 
        r, T, 
        getAssetPrice,
        payoff,
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    
    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myGamma[i]==Approx(BSGamma(myXDomain[i], K, sig, r, T)).epsilon(.0001));
    }
}


TEST_CASE("FSTSPut", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    auto payoff=[&](const auto& assetValue){
        return assetValue<K?K-assetValue:0.0;
    };
    auto getAssetPrice=[&](const auto& logAsset){
        return K*exp(logAsset);
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FSTSPrice(numX, xmax, 
        r, T, 
        getAssetPrice,
        payoff,
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BSPut(myXDomain[i], discount, K, sig, T)).epsilon(.0001));
    }
}

TEST_CASE("FSTSPutDelta", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    auto payoff=[&](const auto& assetValue){
        return assetValue<K?K-assetValue:0.0;
    };
    auto getAssetPrice=[&](const auto& logAsset){
        return K*exp(logAsset);
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myDelta=optionprice::FSTSDelta(numX, xmax, 
        r, T, 
        getAssetPrice,
        payoff,
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myDelta[i]==Approx(BSPutDelta(myXDomain[i], discount, K, sig, T)).epsilon(.0001));
    }
}
TEST_CASE("FSTSPutTheta", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    auto payoff=[&](const auto& assetValue){
        return assetValue<K?K-assetValue:0.0;
    };
    auto getAssetPrice=[&](const auto& logAsset){
        return K*exp(logAsset);
    };

    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myTheta=optionprice::FSTSTheta(numX, xmax, 
        r, T, 
        getAssetPrice,
        payoff,
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    
    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myTheta[i]==Approx(BSPutTheta(myXDomain[i], K, sig, r, T)).epsilon(.0001));
    }
}

TEST_CASE("FSTSPutLowT", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=.25;
    auto K=50.0; 
    
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    auto payoff=[&](const auto& assetValue){
        return assetValue<K?K-assetValue:0.0;
    };
    auto getAssetPrice=[&](const auto& logAsset){
        return K*exp(logAsset);
    };
    int numX=pow(2, 10);
    double xmax=5.0*sqrt(T);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FSTSPrice(numX, xmax, 
        r, T, 
        getAssetPrice,
        payoff,
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BSPut(myXDomain[i], discount, K, sig, T)).epsilon(.0001));
    }
}


TEST_CASE("FangOosterleeCall", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticOption=BSCall(S0, discount, KArray[i], sig, T);
        REQUIRE(myOptionsPrice[i]==Approx(analyticOption).epsilon(.0001));
    }
}
TEST_CASE("FangOosterleeCallDelta", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsDelta=optionprice::FangOostCallDelta(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticDelta=BSCallDelta(S0, discount, KArray[i], sig, T);
        REQUIRE(myOptionsDelta[i]==Approx(analyticDelta).epsilon(.0001));
    }
}
TEST_CASE("FangOosterleeCallTheta", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);
    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsTheta=optionprice::FangOostCallTheta(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticTheta=BSCallTheta(S0, KArray[i], sig, r, T);
        REQUIRE(myOptionsTheta[i]==Approx(analyticTheta).epsilon(.0001));
    }
}
TEST_CASE("FangOosterleeCallGamma", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);
    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsGamma=optionprice::FangOostCallGamma(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticGamma=BSGamma(S0, KArray[i], sig, r, T);
        REQUIRE(myOptionsGamma[i]==Approx(analyticGamma).epsilon(.0001));
    }
}
//rather unecessary, but whatever
TEST_CASE("FangOosterleePutGamma", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);
    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsGamma=optionprice::FangOostPutGamma(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticGamma=BSGamma(S0, KArray[i], sig, r, T);
        REQUIRE(myOptionsGamma[i]==Approx(analyticGamma).epsilon(.0001));
    }
}
TEST_CASE("FangOosterleePutTheta", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsTheta=optionprice::FangOostPutTheta(S0, KArray,  r, T,numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticTheta=BSPutTheta(S0, KArray[i], sig, r, T);
        REQUIRE(myOptionsTheta[i]==Approx(analyticTheta).epsilon(.0001));
    }
}
TEST_CASE("FangOosterleePutDelta", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);
    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsDelta=optionprice::FangOostPutDelta(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticDelta=BSPutDelta(S0, discount, KArray[i], sig, T);
        REQUIRE(myOptionsDelta[i]==Approx(analyticDelta).epsilon(.0001));
    }
}



TEST_CASE("FangOosterleeCallFewK", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);
    int numU=64; //This should run extremely quickly since only have 5 strikes
    auto KArray=std::vector<double>(5, 0);
    KArray[4]=.3;
    KArray[3]=40;//strike 40
    KArray[2]=50;//strike 50
    KArray[1]=60;//strike 60
    KArray[0]=7500;

    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    for(int i=1;i<4; ++i){ 
        auto analyticOption=BSCall(S0, discount, KArray[i], sig, T);
        REQUIRE(myOptionsPrice[i]==Approx(analyticOption).epsilon(.0001));
    }
}
TEST_CASE("FangOosterleeCallFewKLowT", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=.25;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        return exp(chfunctions::gaussLogCF(u, r-sig*sig*.5, sig)*T);
    };
    auto discount=exp(-r*T);
    int numU=256; //This should run extremely quickly since only have 5 strikes, the higher U is to compensate for low t
    auto KArray=std::deque<double>(5, 0);
    KArray[4]=.3;
    KArray[3]=40;//strike 40
    KArray[2]=50;//strike 50
    KArray[1]=60;//strike 60
    KArray[0]=7500;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    for(int i=1;i<4; ++i){ 
        auto analyticOption=BSCall(S0, discount, KArray[i], sig, T);
        REQUIRE(myOptionsPrice[i]==Approx(analyticOption).epsilon(.0001));
    }
}


TEST_CASE("FangOosterleePut", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FangOostPutPrice(S0, KArray, r, T, numU, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticOption=BSPut(S0, discount, KArray[i], sig, T);
        REQUIRE(myOptionsPrice[i]==Approx(analyticOption).epsilon(.0001));
    }
}

TEST_CASE("CarrMadanCall", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0;
    auto discount=exp(-r*T);
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };

    int numX=pow(2, 10);
    auto ada=.25;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::CarrMadanCall(numX,  
        ada,
        S0, 
        r, T, 
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Carr Madan time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    auto b=optionprice::getMaxK(ada);
    auto lambda=optionprice::getLambda(numX, b);
    auto myXDomain=optionprice::getCarrMadanStrikes(ada, S0, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BSCall(S0, discount, myXDomain[i], sig, T)).epsilon(.0001));
    }
}
TEST_CASE("CarrMadanCallLowT", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=.25;
    auto S0=50.0;
    auto discount=exp(-r*T);
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    int numX=pow(2, 10);
    auto ada=.25;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::CarrMadanCall(numX,  
        ada,
        S0, 
        r, T, 
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Carr Madan time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    auto b=optionprice::getMaxK(ada);
    auto lambda=optionprice::getLambda(numX, b);
    auto myXDomain=optionprice::getCarrMadanStrikes(ada, S0, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BSCall(S0, discount, myXDomain[i], sig, T)).epsilon(.0001));
    }
}


TEST_CASE("CarrMadanPut", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0;
    auto discount=exp(-r*T);
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    int numX=pow(2, 10);
    auto ada=.25;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        r, T, 
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Carr Madan time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    auto b=optionprice::getMaxK(ada);
    auto lambda=optionprice::getLambda(numX, b);
    auto myXDomain=optionprice::getCarrMadanStrikes(ada, S0, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BSPut(S0, discount, myXDomain[i], sig, T)).epsilon(.0001));
    }
}



TEST_CASE("FangOosterleeCallCGMY", "[OptionPricing]"){
    //https://cs.uwaterloo.ca/~paforsyt/levy.pdf pg 19
    //S K T r q σ C G M Y
    //90 98 0.25 0.06 0.0 0.0 16.97 7.08 29.97 0.6442
    auto KArray=std::vector<double>(3, 0);
    KArray[2]=.3;
    KArray[1]=98;//strike 98
    KArray[0]=7500;

    auto r=.06;
    auto sig=0.0;
    auto T=.25;
    auto S0=90.0;  
    auto C=16.97;
    auto G=7.08;
    auto M=29.97;
    auto Y=0.6442;

    auto BSCF=[&](const auto& u){
        return exp(chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T);
    };
    
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, BSCF);
    auto myReference=16.212478;//https://cs.uwaterloo.ca/~paforsyt/levy.pdf pg 19
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));


}
TEST_CASE("FangOosterleeCallCGMYWitht=1", "[OptionPricing]"){
    //http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf pg 19
    //S0 = 100, K = 100, r = 0.1, q = 0, C = 1, G = 5, M = 5, T = 1
    //= 19.812948843
    auto KArray=std::vector<double>(3, 0);
    KArray[2]=.3;
    KArray[1]=100.0;//strike 98
    KArray[0]=7500;

    auto r=.1;
    auto sig=0.0;
    auto T=1.0;
    auto S0=100.0;  
    auto C=1.0;
    auto G=5.0;
    auto M=5.0;
    auto Y=0.5;

    auto BSCF=[&](const auto& u){
        return exp(chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T);
    };
    
    int numU=64; //
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, BSCF);
    auto myReference=19.812948843;
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));
}


TEST_CASE("FangOosterleeCallCGMYReducesToBS", "[OptionPricing]"){
    //https://cs.uwaterloo.ca/~paforsyt/levy.pdf pg 20
    //S K T r q σ C G M Y
    //500 500 0.25 0.4 0.0 0.2 1.0 1.4 2.5 1.4
    auto KArray=std::vector<double>(3, 0);
    KArray[2]=.3;
    KArray[1]=500.0;//strike 500.0
    KArray[0]=7500;

    auto r=.4;//seems high
    auto sig=0.2;
    auto T=.25;
    auto S0=500.0;  
    auto C=0.0; //THIS IS ZERO, to reduce to bS
    auto G=1.4;
    auto M=2.5;
    auto Y=1.4;
    auto discount=exp(-r*T);
    auto BSCF=[&](const auto& u){
        return exp(chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T);
    };
    
    int numU=256; //higher due to low T
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, BSCF);
    auto myReference=BSCall(S0, discount, KArray[1], sig, T);
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));
}


TEST_CASE("FangOosterleePutCGMYWithVol", "[OptionPricing]"){
    //https://cs.uwaterloo.ca/~paforsyt/levy.pdf pg 20
    //S K T r q σ C G M Y
    //500 500 0.25 0.4 0.0 0.2 1.0 1.4 2.5 1.4
    auto KArray=std::vector<double>(3, 0);
    KArray[2]=.3;
    KArray[1]=500.0;//strike 500.0
    KArray[0]=7500;

    auto r=.4;//seems high
    auto sig=0.2;
    auto T=.25;
    auto S0=500.0;  
    auto C=1.0;
    auto G=1.4;
    auto M=2.5;
    auto Y=1.4;
    auto discount=exp(-r*T);

    auto BSCF=[&](const auto& u){
        return exp(chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T);
    };
    
    int numU=256; //this is high...this seems to be a computationally tricky problem
    auto myOptionsPrice=optionprice::FangOostPutPrice(S0, KArray, r, T, numU, BSCF);
    //auto myReference=108.49914;//https://cs.uwaterloo.ca/~paforsyt/levy.pdf pg 20
    auto myReference=108.445317144;//This is what I get when I keep increasing U
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));


}

TEST_CASE("CarrMadanCGMY", "[OptionPricing]"){
    //https://cs.uwaterloo.ca/~paforsyt/levy.pdf pg 20
    //S K T r q σ C G M Y
    //500 500 0.25 0.4 0.0 0.2 1.0 1.4 2.5 1.4
    auto r=.4;//seems high
    auto sig=0.2;
    auto T=.25;
    auto S0=500.0;  
    auto C=1.0;
    auto G=1.4;
    auto M=2.5;
    auto Y=1.4;

    auto BSCF=[&](const auto& u){
        return exp(chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T);
    };
    
    
    int numX=pow(2, 10); //this is high...this seems to be a computationally tricky problem
    auto ada=.25;
    auto myOptionsPrice=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        r, T, 
        BSCF
    );
   
    auto b=optionprice::getMaxK(ada);
    auto lambda=optionprice::getLambda(numX, b);
    auto myXDomain=optionprice::getCarrMadanStrikes(ada, S0, numX);
    //auto myReference=108.49914;//https://cs.uwaterloo.ca/~paforsyt/levy.pdf pg 20
    auto myReference=108.4614527579;//This is what I get when I keep increasing U 
    REQUIRE(myOptionsPrice[numX/2]==Approx(myReference).epsilon(.0001));

}




//correlation on diffusion only

TEST_CASE("CarrMadanCGMYPut", "[OptionPricing]"){
    auto r=.04;//seems high
    auto sig=0.2;
    auto T=.25;
    auto S0=500.0;  
    auto C=1.0;
    auto G=1.4;
    auto M=2.5;
    auto Y=.6;
    auto kappa=.5;
    auto a=kappa;//.5;//long run tau of 1
    auto v0=1.05;
    auto adaV=.2;
    auto rho=-1.0;//-.4;//leverage rho
    auto CFCorr=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::cgmyLogRNCF(u, C, G, M, Y, 0.0, sig),
            a, 
            kappa-adaV*rho*u*sig,
            adaV,
            T, 
            v0
        ));
    };
    auto CFCorrOnlyBM=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::gaussLogCF(u, -sig*sig*.5, sig),
            a, 
            kappa-adaV*rho*u*sig,
            adaV,
            T, 
            v0
        )+chfunctions::cgmyLogRNCF(u, C, G, M, Y, 0.0, 0.0)*T);
    };
    auto CFStoch=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::cgmyLogRNCF(u, C, G, M, Y, 0.0, sig),
            a, 
            kappa,
            adaV,
            T, 
            v0
        ));
    };
    auto CFBase=[&](const auto& u){
        return exp(
            chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T
        );
    };
    int numX=pow(2, 10); 
    auto ada=.25;
    auto myOptionsPriceCorr=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        r, T, 
        CFCorr
    );
    auto myOptionsPriceStoch=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        r, T, 
        CFStoch
    );
    auto myOptionsPriceBase=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        r, T, 
        CFBase
    );
    auto myOptionsPriceCorrBM=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        r, T, 
        CFCorrOnlyBM
    );
   
    auto b=optionprice::getMaxK(ada);
    auto lambda=optionprice::getLambda(numX, b);
    auto myXDomain=optionprice::getCarrMadanStrikes(ada, S0, numX);
    auto myReferenceCorr=74.9251;//I think this is right, but this is based solely off whether it makes logical sense
    auto myReferenceStoch=75.0101;//I think this is right, but this is based solely off whether it makes logical sense
    auto myReferenceBase=72.9703833045;//I think this is right, but this is based solely off whether it makes logical sense
    std::cout<<"corr: "<<myOptionsPriceCorr[numX/2]<<std::endl;
    std::cout<<"stoch: "<<myOptionsPriceStoch[numX/2]<<std::endl;
    std::cout<<"base: "<<myOptionsPriceBase[numX/2]<<std::endl;
    std::cout<<"corrBM: "<<myOptionsPriceCorrBM[numX/2]<<std::endl;
    std::cout<<"corr: "<<myOptionsPriceCorr[numX/2-100]<<std::endl;
    std::cout<<"stoch: "<<myOptionsPriceStoch[numX/2-100]<<std::endl;
    std::cout<<"base: "<<myOptionsPriceBase[numX/2-100]<<std::endl;
    std::cout<<"corrBM: "<<myOptionsPriceCorrBM[numX/2-100]<<std::endl;
    REQUIRE(myOptionsPriceCorr[numX/2]==Approx(myReferenceCorr).epsilon(.0001));
    REQUIRE(myOptionsPriceStoch[numX/2]==Approx(myReferenceStoch).epsilon(.0001));
    REQUIRE(myOptionsPriceBase[numX/2]==Approx(myReferenceBase).epsilon(.0001));

}

/**test heston*/
//http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf pg 15
TEST_CASE("CarrMadanHeston", "[OptionPricing]"){
    //heston parameters
    double b=.0398;
    double a=1.5768;
    double c=.5751;
    double rho=-.5711;
    double v0=.0175;
    auto r=0.0;
    //convert to extended CGMY
    auto sig=sqrt(b);
    auto speed=a;
    auto T=1.0;
    auto S0=100.0;  
    auto kappa=speed;//long run tau of 1
    auto v0Hat=v0/b;
    auto adaV=c/sqrt(b);

    /**These two CFs should be the same*/
    auto CFCorr=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::cgmyLogRNCF(u, 0.0, 1.0, 1.0, .5, 0.0, sig),
            speed, 
            kappa-adaV*rho*u*sig,
            adaV,
            T, 
            v0Hat
        ));
    };
    auto CFCorrOnlyBM=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::gaussLogCF(u, -sig*sig*.5, sig),
            speed, 
            kappa-adaV*rho*u*sig,
            adaV,
            T, 
            v0Hat
        ));
    };
    int numX=pow(2, 10); 
    auto ada=.25;
    auto myOptionsPriceCorr=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        r, T, 
        CFCorr
    );

    auto myOptionsPriceCorrBM=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        r, T, 
        CFCorrOnlyBM
    );
   
    auto bLim=optionprice::getMaxK(ada);
    auto lambda=optionprice::getLambda(numX, bLim);
    auto myXDomain=optionprice::getCarrMadanStrikes(ada, S0, numX);
    auto myReference= 5.78515545;
    std::cout<<myOptionsPriceCorr[numX/2]<<std::endl;
    std::cout<<myOptionsPriceCorrBM[numX/2]<<std::endl;
    REQUIRE(myOptionsPriceCorr[numX/2]==Approx(myReference).epsilon(.0001));
    REQUIRE(myOptionsPriceCorrBM[numX/2]==Approx(myReference).epsilon(.0001));
}
TEST_CASE("FangOostHeston", "[OptionPricing]"){
    //heston parameters
    double b=.0398;
    double a=1.5768;
    double c=.5751;
    double rho=-.5711;
    double v0=.0175;
    auto r=0.0;
    //convert to extended CGMY
    auto sig=sqrt(b);
    auto speed=a;
    auto T=1.0;
    auto S0=100.0;  
    auto kappa=speed;//long run tau of 1
    auto v0Hat=v0/b;
    auto adaV=c/sqrt(b);
    auto CFCorr=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::cgmyLogRNCF(u, 0.0, 1.0, 1.0, .5, 0.0, sig),
            speed, 
            kappa-adaV*rho*u*sig,
            adaV,
            T, 
            v0Hat
        ));
    };
    std::vector<double> KArray(3);
    KArray[2]=.3;
    KArray[1]=100;
    KArray[0]=5000;
    int numU=256; //this is high...this seems to be a computationally tricky problem
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, CFCorr);
    auto myReference= 5.78515545;
    std::cout<<myOptionsPrice[1]<<std::endl;
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));
}
TEST_CASE("FangOostHestonSubsetMerton", "[OptionPricing]"){
    /**note that the only difference between this and the other Heston calculations is
    that I'm using Merton's CF*/
    //heston parameters
    double b=.0398;
    double a=1.5768;
    double c=.5751;
    double rho=-.5711;
    double v0=.0175;
    auto r=0.0;
    //convert to extended CGMY
    auto sig=sqrt(b);
    auto speed=a;
    auto T=1.0;
    auto S0=100.0;  
    auto kappa=speed;//long run tau of 1
    auto v0Hat=v0/b;
    auto adaV=c/sqrt(b);
    /**const T& u, const Number& lambda, const Number& muL, const Number& sigL, const Number& r,  const Number& sigma*/
    auto CFCorr=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::mertonLogRNCF(u, 0.0, 1.0, 1.0, 0.0, sig),
            speed, 
            kappa-adaV*rho*u*sig,
            adaV,
            T, 
            v0Hat
        ));
    };
    std::vector<double> KArray(3);
    KArray[2]=.3;
    KArray[1]=100;
    KArray[0]=5000;
    int numU=256; //this is high...this seems to be a computationally tricky problem
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, CFCorr);
    auto myReference= 5.78515545;
    std::cout<<myOptionsPrice[1]<<std::endl;
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));
}
//https://www.upo.es/personal/jfernav/papers/Jumps_JOD_.pdf pg 8
TEST_CASE("FangOostMerton", "[OptionPricing]"){
    //merton parameters
    double r=.1;
    double lambda=1;
    double sigmaL=sqrt(.05);
    double sigma=sigmaL;
    double mu=-sigmaL*sigmaL*.5;
    double S0=38;
    double T=.5;
    auto CFCorr=[&](const auto& u){
        return exp(chfunctions::mertonLogRNCF(u, lambda, mu, sigmaL, r, sigma)*T);
    };
    std::vector<double> KArray(3);
    KArray[2]=.3;
    KArray[1]=35;
    KArray[0]=5000;
    int numU=256; //
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, CFCorr);
    auto myReference= 5.9713;
    std::cout<<myOptionsPrice[1]<<std::endl;
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));
}

TEST_CASE("FSTSHeston", "[OptionPricing]"){
    //heston parameters
    double b=.0398;
    double a=1.5768;
    double c=.5751;
    double rho=-.5711;
    double v0=.0175;
    auto r=0.0;
    //convert to extended CGMY
    auto sig=sqrt(b);
    auto speed=a;
    auto T=1.0;
    //auto S0=100.0;  
    auto strike=100.0;
    auto kappa=speed;//long run tau of 1
    auto v0Hat=v0/b;
    auto adaV=c/sqrt(b);
    auto getAssetPrice=[&](const auto& logAsset){
        return strike*exp(logAsset);
    };
    auto CFCorr=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::cgmyLogRNCF(u, 0.0, 1.0, 1.0, .5, 0.0, sig),
            speed, 
            kappa-adaV*rho*u*sig,
            adaV,
            T, 
            v0Hat
        ));
    };
    int numU=256; //this is high...this seems to be a computationally tricky problem
    auto payoff=[&](const auto& assetValue){
        return assetValue>strike?assetValue-strike:0.0;
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FSTSPrice(numX, xmax, 
        r, T,
        getAssetPrice, 
        payoff, 
        CFCorr
    );
    auto myXDomain=optionprice::getStrikeUnderlying(-xmax, xmax, strike, numX);
    std::cout<<myXDomain[numX/2]<<std::endl;
    auto myReference= 6.09466;
    std::cout<<myOptionsPrice[numX/2]<<std::endl;
    REQUIRE(myOptionsPrice[numX/2]==Approx(myReference).epsilon(.001));
}



TEST_CASE("FangOostDegenerateBS", "[OptionPricing]"){
    auto sig=.3;
    auto r=.03;
    auto adaV=0.0;//zero, so should be BS
    auto speed=.4;//shouldn't matter
    auto T=.25;
    auto rho=.3;//shouldn't matter
    auto discount=exp(-r*T);
    auto S0=100.0;
    auto v0Hat=1.0;//else will be weird..a decay of a non-stochastic time change to one
    auto CFCorr=[&](const auto& u){
        return exp(r*T*u+chfunctions::cirLogMGF(
            -chfunctions::cgmyLogRNCF(u, 0.0, 1.0, 1.0, .5, 0.0, sig),
            speed, 
            speed-adaV*rho*u*sig,
            adaV,
            T, 
            v0Hat
        ));
    };
    std::vector<double> KArray(3); 
    KArray[2]=.3;
    KArray[1]=100;
    KArray[0]=5000;
    int numU=256;
    auto myOptionsPrice=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, CFCorr);
    auto myReference= BSCall(S0, discount, KArray[1], sig, T);
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));
}

TEST_CASE("transformPrices", "[OptionCalibration]"){
    std::vector<double> arr={3, 4, 5};
    double asset=4.0;
    double minV=2.0;
    double maxV=6.0;
    auto result=optioncal::transformPrices(arr, asset, minV, maxV);
    std::vector<double> expected={.5, .75, 1.0, 1.25, 1.5};
    REQUIRE(result==expected);
}

TEST_CASE("filter returns empty first array when condition never met", "[OptionCalibration]"){
    std::vector<double> arr={3, 4, 5};
    std::vector<double> expected1={};
    std::vector<double> expected2=arr;
    auto result=optioncal::filter(arr, [](const auto& x, const auto& index){
            return false;
        }, 
        [](const auto& x, const auto& index){return x;}, 
        [](const auto& x, const auto& index){return x;}
    );
    REQUIRE(std::get<0>(result)==expected1);
    REQUIRE(std::get<1>(result)==expected2);
}

TEST_CASE("filter returns two arrays when condition sometimes met", "[OptionCalibration]"){
    std::vector<double> arr={3, 4, 5};
    std::vector<double> expected1={3, 5};
    std::vector<double> expected2={4};
    auto result=optioncal::filter(arr, [](const auto& x, const auto& index){
        return index%2==0;
    }, [](const auto& x, const auto& index){return x;}, [](const auto& x, const auto& index){return x;});
    REQUIRE(std::get<0>(result)==expected1);
    REQUIRE(std::get<1>(result)==expected2);
}
TEST_CASE("getObjfn_arr", "[OptionCalibration]"){
    std::vector<std::complex<double> > arr={
        std::complex<double>(3, 0), 
        std::complex<double>(4, 0), 
        std::complex<double>(5, 0)
    };
    std::vector<double> u={6, 7, 8};
    auto cf=[](const auto& _){
        return [](const auto& cmpU){
            return std::complex<double>(cmpU.imag(), 0.0);
        };
    };
    auto hoc=optioncal::getObjFn_arr(std::move(arr), std::move(cf), std::move(u));
    double expected=9;//3*3^2/3
    double tmp=0;
    REQUIRE(hoc(tmp)==expected);

}

TEST_CASE("getObjfn_arr multiple parameters", "[OptionCalibration]"){
    std::vector<std::complex<double> > arr={
        std::complex<double>(3, 0), 
        std::complex<double>(4, 0), 
        std::complex<double>(5, 0)
    };
    std::vector<double> u={6, 7, 8};
    auto cf=[](const auto& v1, const auto& v2){
        return [&](const auto& cmpU){
            return std::complex<double>(cmpU.imag(), 0.0);
        };
    };
    /*auto cf=[](const auto& cmpU, const auto& v1, const auto& v2){
        return std::complex<double>(cmpU.imag(), 0.0);
    };*/
    auto hoc=optioncal::getObjFn_arr(std::move(arr), std::move(cf), std::move(u));
    double expected=9;//3*3^2/3
    double tmp1=0;
    double tmp2=0;
    REQUIRE(hoc(tmp1, tmp2)==expected);
}
TEST_CASE("getObjfn_arr vector parameter", "[OptionCalibration]"){
    std::vector<std::complex<double> > arr={
        std::complex<double>(3, 0), 
        std::complex<double>(4, 0), 
        std::complex<double>(5, 0)
    };
    std::vector<double> u={6, 7, 8};
    /*auto cf=[](const auto& cmpU, const auto& v1){
        return std::complex<double>(cmpU.imag(), 0.0);
    };*/
    auto cf=[](const auto& v1){
        return [&](const auto& cmpU){
            return std::complex<double>(cmpU.imag(), 0.0);
        };
    };
    auto hoc=optioncal::getObjFn_arr(std::move(arr), std::move(cf), std::move(u));
    double expected=9;//3*3^2/3
    std::vector<double> tmp1={3, 4, 5};
    REQUIRE(hoc(tmp1)==expected);
}
TEST_CASE("getKThatIsBelowOne with exactly 1", "[OptionCalibration]"){
    std::vector<double> arr={.8, .9, 1, 1.1, 1.2};
    int expected=2;
    REQUIRE(optioncal::getKThatIsBelowOne(arr)==expected);
}
TEST_CASE("getKThatIsBelowOne with no exact 1", "[OptionCalibration]"){
    std::vector<double> arr={.8, .9, 1.1, 1.2};
    int expected=1;
    REQUIRE(optioncal::getKThatIsBelowOne(arr)==expected);
}