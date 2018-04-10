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

/*
TEST_CASE("test xJ", "[OptionCalibration]"){
    double maxStrike=100;
    double stock=10;
    std::cout<<"min xJ: "<<optioncal::xJ(stock, stock/maxStrike, 1.0)<<std::endl;
    std::cout<<"max xJ: "<<optioncal::xJ(stock, stock*maxStrike, 1.0)<<std::endl;
    
}
TEST_CASE("test fSpline", "[OptionCalibration]"){
    double stock=10.0;
    double discount=.99;
    double sigma=.3;
    double T=1.0;
    double r=-log(discount)/T;
    std::vector<double> strikes={4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    auto optionPrices=futilities::for_each_parallel_copy(strikes, [&](const auto& k, const auto& index){
        return BSCall(stock, discount, k, sigma, T);
    });
    double maxStrike=100;
    auto estimateOfPhi=optioncal::generateFOEstimate(strikes, optionPrices, stock, discount, maxStrike);
    auto result=estimateOfPhi(1.0);
    std::cout<<"result of phi at 1.0: "<<result<<std::endl;
}*/

TEST_CASE("objFn", "[OptionCalibration]"){
    double stock=10.0;
    double discount=.99;
    double sigma=.3;
    double T=1.0;
    double r=-log(discount)/T;
    std::vector<double> strikes={4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    auto optionPrices=futilities::for_each_parallel_copy(strikes, [&](const auto& k, const auto& index){
        return BSCall(stock, discount, k, sigma, T);
    });
    double maxStrike=100;
    auto estimateOfPhi=optioncal::generateFOEstimate(strikes, optionPrices, stock, discount, maxStrike);
    //int N=20;
    
    //auto myCf=
    auto getU=[](const auto& N){
        double du= 2*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    auto objFn=optioncal::getObjFn(
        std::move(estimateOfPhi),
        [&](const auto& u, const auto& sigma){
            return u*(r-sigma*sigma*.5)*T+sigma*sigma*.5*T*u*u;
        },
        getU(30)
    );
    std::cout<<"guess=.4: "<<objFn(.4)<<std::endl;
    std::cout<<"guess=.3: "<<objFn(.3)<<std::endl;
}



TEST_CASE("BSCal", "[OptionCalibration]"){
    double stock=10.0;
    double discount=1.0; //this really impacts the accuracy
    double sigma=.3;
    double T=1.0;
    double r=-log(discount)/T;
    std::vector<double> strikes={4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    auto optionPrices=futilities::for_each_parallel_copy(strikes, [&](const auto& k, const auto& index){
        return BSCall(stock, discount, k, sigma, T);
    });
    double minStrike=.025;
    double maxStrike=200;
    auto estimateOfPhi=optioncal::generateFOEstimate(strikes, optionPrices, stock, discount, maxStrike);
    auto estimateOfPhiSpline=optioncal::generateFOEstimateSpline(strikes, optionPrices, stock, discount, minStrike, maxStrike);
    //int N=20;
    
    //auto myCf=
    /*auto getU=[](const auto& N){
        double maxU=20.0;
        double du= 2*maxU/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return -maxU+index*du;
        });
    };*/
    /**theoretically I should have this on -infty to infty...but it seems to work on 0 to 2pi*/
    auto getU=[](const auto& N){
        double du= 2*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    
    auto tmpCF=[&](const auto& u){
        return u*(r-sigma*sigma*.5)*T+sigma*sigma*.5*T*u*u;
    };

    /*std::cout<<"BS"<<std::endl;
    std::cout<<"phiHat @ -100: "<<estimateOfPhi(-100.0)<<", cf @ -100: "<<tmpCF(std::complex<double>(0.0, -100.0))<<std::endl;
    std::cout<<"phiHat @ -2: "<<estimateOfPhi(-2.0)<<", cf @ -2: "<<tmpCF(std::complex<double>(0.0, -2.0))<<std::endl;
    std::cout<<"phiHat @ -.5: "<<estimateOfPhi(-0.5)<<", cf @ -.5: "<<tmpCF(std::complex<double>(0.0, -.5))<<std::endl;
    std::cout<<"phiHat @ 0: "<<estimateOfPhi(0.0)<<", cf @ 0: "<<tmpCF(std::complex<double>(0.0, 0.0))<<std::endl;
    std::cout<<"phiHat @ .5: "<<estimateOfPhi(.5)<<", cf @ .5: "<<tmpCF(std::complex<double>(0.0, .5))<<std::endl;
    std::cout<<"phiHat @ 2: "<<estimateOfPhi(2.0)<<", cf @ 2: "<<tmpCF(std::complex<double>(0.0, 2.0))<<std::endl;
    std::cout<<"phiHat @ 100: "<<estimateOfPhi(100.0)<<", cf @ 100: "<<tmpCF(std::complex<double>(0.0, 100.0))<<std::endl;*/

    
    /*std::cout<<"BS with FFT"<<std::endl;
    int N=128;
    double xMin=log(minStrike*discount/stock);
    double xMax=log(maxStrike*discount/stock);
    auto phis=estimateOfPhiSpline(N);
    double uMax=M_PI*(N-1)/(xMax-xMin);
    double du=2.0*uMax/N;
    for(int i=0; i<phis.size(); ++i){ //*std::complex<double>(0, 1.0)+1.0
        std::cout<<"u: "<<i*du<<", estimate: "<<phis[i]<<", exact: "<<tmpCF(std::complex<double>(0.0, i*du))<<std::endl;
    }
*/
   
    std::cout<<"BS with FFT"<<std::endl;
    const auto estimateOfPhiDHS=optioncal::generateFOEstimateDHS(strikes, optionPrices, discount, stock, minStrike, maxStrike);
    int N=128;
    double xMin=log(minStrike/stock);
    double xMax=log(maxStrike/stock);
    auto phis=estimateOfPhiDHS(N);
    //double uMin=-20.0;
    double uMax=M_PI*(N-1)/(xMax-xMin);
    double du=(2.0*uMax)/N; 
    for(int i=0; i<phis.size(); ++i){ //*std::complex<double>(0, 1.0)+1.0
        std::cout<<"u: "<<i*du-uMax<<", estimate: "<<phis[i]<<", exact: "<<tmpCF(std::complex<double>(0.0, i*du-uMax))<<std::endl;
    }

    auto objFn=optioncal::getObjFn(
        std::move(estimateOfPhi),
        [&](const auto& u, const auto& sigma){
            return u*(r-sigma*sigma*.5)*T+sigma*sigma*.5*T*u*u;
        },
        getU(30)
    ); 
    auto guess=.2;
    auto results=optioncal::calibrate(objFn, guess);
    std::cout<<"optimal sigma: "<<std::get<0>(results)<<std::endl;
    std::cout<<"value at actual: "<<objFn(sigma)<<std::endl;
    
}



TEST_CASE("HestonCal", "[OptionCalibration]"){
    
    double b=.0398;
    double a=1.5768;
    double c=.5751;
    double rho=-.5711;
    double v0=.0175;
    auto r=0.01;
    //convert to extended CGMY
    auto sig=sqrt(b);
    auto speed=a;
    auto T=1.0;
    auto S0=100.0;  
    //auto kappa=speed;//long run tau of 1
    auto v0Hat=v0/b;
    auto adaV=c/sqrt(b);
    auto CFCorr=[&](
        const auto& u,
        const auto& sig, 
        const auto& speed, 
        const auto& adaV,
        const auto& rho, 
        const auto& v0Hat
    ){
        return r*T*u+chfunctions::cirLogMGF(
            -chfunctions::cgmyLogRNCF(u, 0.0, 1.0, 1.0, .5, 0.0, sig),
            speed, 
            speed-adaV*rho*u*sig,
            adaV,
            T, 
            v0Hat
        );
    };
    std::vector<double> KArray={
        5000, 130, 120, 110, 105, 100, 95, 90, 80, 70, .3
    };
    int numU=256; //this is high...this seems to be a computationally tricky problem
    auto optionPrices=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, [&](const auto& u){
        return exp(CFCorr(u, sig, speed, adaV, rho, v0Hat));
    });


    double discount=exp(-r*T);
    double maxStrike=160;
    auto getReducedAndReversed=[](const std::vector<double>& arr){
        auto arrTmp=std::vector<double>(arr.begin()+1, arr.end()-1);
        std::reverse(std::begin(arrTmp), std::end(arrTmp));
        return arrTmp;
    };
    auto observedK=getReducedAndReversed(KArray);
    auto observedO=getReducedAndReversed(optionPrices);

    for(auto& k:observedK){
        std::cout<<"k: "<<k<<std::endl;
    }
    for(auto& c:observedO){
        std::cout<<"c: "<<c<<std::endl;
    }
    
    auto estimateOfPhi=optioncal::generateFOEstimate(
        observedK, 
        observedO, 
        S0, discount, maxStrike
    );

    auto getU=[](const auto& N){
        double du= 2*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    std::cout<<"Heston "<<std::endl;
    std::cout<<"phiHat @ -100: "<<estimateOfPhi(-100.0)<<", cf @ -100: "<<CFCorr(std::complex<double>(0.0, -100.0), sig, speed, adaV, rho, v0Hat)<<std::endl;
    std::cout<<"phiHat @ -2: "<<estimateOfPhi(-2.0)<<", cf @ -2: "<<CFCorr(std::complex<double>(0.0, -2.0), sig, speed, adaV, rho, v0Hat)<<std::endl;
    std::cout<<"phiHat @ -.5: "<<estimateOfPhi(-0.5)<<", cf @ -.5: "<<CFCorr(std::complex<double>(0.0, -.5), sig, speed, adaV, rho, v0Hat)<<std::endl;
    std::cout<<"phiHat @ 0: "<<estimateOfPhi(0.0)<<", cf @ 0: "<<CFCorr(std::complex<double>(0.0, 0.0), sig, speed, adaV, rho, v0Hat)<<std::endl;
    std::cout<<"phiHat @ .5: "<<estimateOfPhi(.5)<<", cf @ .5: "<<CFCorr(std::complex<double>(0.0, .5), sig, speed, adaV, rho, v0Hat)<<std::endl;
    std::cout<<"phiHat @ 2: "<<estimateOfPhi(2.0)<<", cf @ 2: "<<CFCorr(std::complex<double>(0.0, 2.0), sig, speed, adaV, rho, v0Hat)<<std::endl;
    std::cout<<"phiHat @ 100: "<<estimateOfPhi(100.0)<<", cf @ 100: "<<CFCorr(std::complex<double>(0.0, 100.0), sig, speed, adaV, rho, v0Hat)<<std::endl;



    auto objFn=optioncal::getObjFn(
        std::move(estimateOfPhi),
        std::move(CFCorr),
        getU(20)
    );
    auto guessSigma=.5; 
    auto guessSpeed=.8;
    auto guessAdaV=1.0; 
    auto guessRho=-.7;
    auto guessV0=.8;



    /*auto results=optioncal::calibrate(objFn, guessSigma, guessSpeed, guessAdaV, guessRho, guessV0);
    std::cout<<"optimal sigma: "<<std::get<0>(results)<<std::endl;
    std::cout<<"optimal speed: "<<std::get<1>(results)<<std::endl;
    std::cout<<"optimal adaV: "<<std::get<2>(results)<<std::endl;
    std::cout<<"optimal rho: "<<std::get<3>(results)<<std::endl;
    std::cout<<"optimal V0: "<<std::get<4>(results)<<std::endl;

    std::cout<<"actual sigma: "<<sig<<std::endl;
    std::cout<<"actual speed: "<<speed<<std::endl;
    std::cout<<"actual adaV: "<<adaV<<std::endl;
    std::cout<<"actual rho: "<<rho<<std::endl;
    std::cout<<"actual v0: "<<v0Hat<<std::endl;
    std::cout<<"obj at optimal: "<<objFn(std::get<0>(results), std::get<1>(results), std::get<2>(results), std::get<3>(results), std::get<4>(results))<<std::endl;

    std::cout<<"obj at actual: "<<objFn(sig, speed, adaV, rho, v0Hat)<<std::endl;*/


}



TEST_CASE("CGMYCal", "[OptionCalibration]"){
    
    auto r=.04;//seems high
    auto sig=0.2;
    auto T=1.0;
    auto S0=100.0;  
    auto C=1.0;
    auto G=1.4;
    auto M=2.5;
    auto Y=.6;

    auto CFBase=[&](const auto& u, const auto& sig, const auto& C, const auto& G, const auto& M, const auto& Y){
        return chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T;
    };
    std::vector<double> KArray={
        5000, 130, 120, 110, 105, 100, 95, 90, 80, 70, .3
    };
    int numU=256; //this is high...this seems to be a computationally tricky problem
    auto optionPrices=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, [&](const auto& u){
        return exp(CFBase(u, sig, C, G, M, Y));
    }); 


    double discount=exp(-r*T);
    double maxStrike=5000;
    auto getReducedAndReversed=[](const std::vector<double>& arr){
        auto arrTmp=std::vector<double>(arr.begin()+1, arr.end()-1);
        std::reverse(std::begin(arrTmp), std::end(arrTmp));
        return arrTmp;
    };
    auto observedK=getReducedAndReversed(KArray);
    auto observedO=getReducedAndReversed(optionPrices);

    for(auto& k:observedK){
        std::cout<<"k: "<<k<<std::endl;
    }
    for(auto& c:observedO){
        std::cout<<"c: "<<c<<std::endl;
    }
    
    auto estimateOfPhi=optioncal::generateFOEstimate(
        observedK, 
        observedO, 
        S0, discount, maxStrike
    );

    auto getU=[](const auto& N){
        double du= 2*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };

    std::cout<<"CGMY"<<std::endl;
    std::cout<<"phiHat @ -100: "<<estimateOfPhi(-100.0)<<", cf @ -100: "<<CFBase(std::complex<double>(0.0, -100.0), sig, C, G, M, Y)<<std::endl;
    std::cout<<"phiHat @ -2: "<<estimateOfPhi(-2.0)<<", cf @ -2: "<<CFBase(std::complex<double>(0.0, -2.0), sig, C, G, M, Y)<<std::endl;
    std::cout<<"phiHat @ -.5: "<<estimateOfPhi(-0.5)<<", cf @ -.5: "<<CFBase(std::complex<double>(0.0, -.5), sig, C, G, M, Y)<<std::endl;
    std::cout<<"phiHat @ 0: "<<estimateOfPhi(0.0)<<", cf @ 0: "<<CFBase(std::complex<double>(0.0, 0.0), sig, C, G, M, Y)<<std::endl;
    std::cout<<"phiHat @ .5: "<<estimateOfPhi(.5)<<", cf @ .5: "<<CFBase(std::complex<double>(0.0, .5), sig, C, G, M, Y)<<std::endl;
    std::cout<<"phiHat @ 2: "<<estimateOfPhi(2.0)<<", cf @ 2: "<<CFBase(std::complex<double>(0.0, 2.0), sig, C, G, M, Y)<<std::endl;
    std::cout<<"phiHat @ 100: "<<estimateOfPhi(100.0)<<", cf @ 100: "<<CFBase(std::complex<double>(0.0, 100.0), sig, C, G, M, Y)<<std::endl;


    auto objFn=optioncal::getObjFn(
        std::move(estimateOfPhi),
        std::move(CFBase),
        getU(20)
    );
    auto guessSigma=.5; 
    auto guessC=.5;
    auto guessG=1.35;
    auto guessM=1.2;
    auto guessY=1.2;

    
    const int maxIter=500;
    const double prec=.00001; 
    const double peterb=.000001;
    const double alpha=.01;//*value; //needs a very small step or it goes off to no where
    /*auto results=newton::gradientDescentApprox([&](const auto& sig, const auto& C, const auto& G, const auto& M, const auto& Y){
        return futilities::sum(optionprice::FangOostCallPrice(S0, KArray, r, T, numU, [&](const auto& u){
            return exp(CFBase(u, sig, C, G, M, Y));
        }), [&](const auto& price, const auto& i){
            return i>0&&i<(optionPrices.size()-1)?futilities::const_power((optionPrices[i]-price)/optionPrices[i], 2):0.0;
        });
    }, maxIter, prec, peterb, alpha, guessSigma, guessC, guessG, guessM, guessY);*/
    



    /*int result=optioncal::calibrate_de([&](const auto& vec){
        return futilities::sum(optionprice::FangOostCallPrice(S0, KArray, r, T, numU, [&](const auto& u){
            return exp(CFBase(u, vec(0), vec(1), vec(2), vec(3), vec(4)));
        }), [&](const auto& price, const auto& i){
            return i>0&&i<(optionPrices.size()-1)?futilities::const_power(optionPrices[i]-price, 2)/optionPrices[i]:0.0;
        });
    }, std::vector<double>({.5, .5, 1.35, 1.2, 1.2}));*/
    
    //int result=optioncal::calibrate_de();




    //auto results=optioncal::calibrate(objFn, guessSigma, guessC, guessG, guessM, guessY);
    /*std::cout<<"optimal sigma: "<<std::get<0>(results)<<std::endl;
    std::cout<<"optimal C: "<<std::get<1>(results)<<std::endl;
    std::cout<<"optimal G: "<<std::get<2>(results)<<std::endl;
    std::cout<<"optimal M: "<<std::get<3>(results)<<std::endl;
    std::cout<<"optimal Y: "<<std::get<4>(results)<<std::endl;

    std::cout<<"actual sigma: "<<sig<<std::endl;
    std::cout<<"actual C: "<<C<<std::endl;
    std::cout<<"actual G: "<<G<<std::endl;
    std::cout<<"actual M: "<<M<<std::endl;
    std::cout<<"actual Y: "<<Y<<std::endl;
    std::cout<<"obj at optimal: "<<objFn(std::get<0>(results), std::get<1>(results), std::get<2>(results), std::get<3>(results), std::get<4>(results))<<std::endl;

    std::cout<<"obj at actual: "<<objFn(sig, C, G, M, Y)<<std::endl;
*/

}

TEST_CASE("Jump diffusion cal", "[OptionCalibration]"){
    
    auto r=.04;//seems high
    auto sig=0.2;
    auto T=.25;
    auto S0=100.0;  
    auto lambda=.05;
    auto muJ=.05;
    auto sigJ=.05;

    auto CFBase=[&](const auto& u, const auto& sig, const auto& lambda, const auto& muJ, const auto& sigJ){
        return (r-sig*sig*.5-lambda*(exp(muJ+sigJ*sigJ*.5)-1.0))*T*u+sig*sig*.5*u*u*T+lambda*(exp(muJ*u+sigJ*sigJ*.5*u*u)-1.0)*T;
    };
    std::vector<double> KArray={
        5000, 130, 120, 110, 105, 100, 95, 90, 80, 70, .3
    };
    int numU=256; //this is high...this seems to be a computationally tricky problem
    auto optionPrices=optionprice::FangOostCallPrice(S0, KArray, r, T, numU, [&](const auto& u){
        return exp(CFBase(u, sig, lambda, muJ, sigJ));
    });


    double discount=exp(-r*T);
    double maxStrike=160;
    auto getReducedAndReversed=[](const std::vector<double>& arr){
        auto arrTmp=std::vector<double>(arr.begin()+1, arr.end()-1);
        std::reverse(std::begin(arrTmp), std::end(arrTmp));
        return arrTmp;
    };
    auto observedK=getReducedAndReversed(KArray);
    auto observedO=getReducedAndReversed(optionPrices);

    for(auto& k:observedK){
        std::cout<<"k: "<<k<<std::endl;
    }
    for(auto& c:observedO){
        std::cout<<"c: "<<c<<std::endl;
    }
    
    auto estimateOfPhi=optioncal::generateFOEstimate(
        observedK, 
        observedO, 
        S0, discount, maxStrike
    );

    auto getU=[](const auto& N){
        double du= 2*M_PI/N;
        return futilities::for_each_parallel(1, N, [&](const auto& index){
            return index*du;
        });
    };
    
    std::cout<<"Jump Diffusion"<<std::endl;

    std::cout<<"phiHat @ -100: "<<estimateOfPhi(-100.0)<<", cf @ -100: "<<CFBase(std::complex<double>(0.0, -100.0), sig, lambda, muJ, sigJ)<<std::endl;
    std::cout<<"phiHat @ -2: "<<estimateOfPhi(-2.0)<<", cf @ -2: "<<CFBase(std::complex<double>(0.0, -2.0), sig, lambda, muJ, sigJ)<<std::endl;
    std::cout<<"phiHat @ -.5: "<<estimateOfPhi(-0.5)<<", cf @ -.5: "<<CFBase(std::complex<double>(0.0, -.5), sig, lambda, muJ, sigJ)<<std::endl;
    std::cout<<"phiHat @ 0: "<<estimateOfPhi(0.0)<<", cf @ 0: "<<CFBase(std::complex<double>(0.0, 0.0), sig, lambda, muJ, sigJ)<<std::endl;
    std::cout<<"phiHat @ .5: "<<estimateOfPhi(.5)<<", cf @ .5: "<<CFBase(std::complex<double>(0.0, .5), sig, lambda, muJ, sigJ)<<std::endl;
    std::cout<<"phiHat @ 2: "<<estimateOfPhi(2.0)<<", cf @ 2: "<<CFBase(std::complex<double>(0.0, 2.0), sig, lambda, muJ, sigJ)<<std::endl;
    std::cout<<"phiHat @ 100: "<<estimateOfPhi(100.0)<<", cf @ 100: "<<CFBase(std::complex<double>(0.0, 100.0), sig, lambda, muJ, sigJ)<<std::endl;

    auto objFn=optioncal::getObjFn(
        std::move(estimateOfPhi),
        std::move(CFBase),
        getU(20)
    );
    auto guessSigma=.3; 
    auto guessLambda=.5;
    auto guessMuJ=.5;
    auto guessSigJ=.5;
   

/*
    auto results=optioncal::calibrate(objFn, guessSigma, guessLambda, guessMuJ, guessSigJ);
    std::cout<<"optimal sigma: "<<std::get<0>(results)<<std::endl;
    std::cout<<"optimal Lambda: "<<std::get<1>(results)<<std::endl;
    std::cout<<"optimal muJ: "<<std::get<2>(results)<<std::endl;
    std::cout<<"optimal sigJ: "<<std::get<3>(results)<<std::endl;

    std::cout<<"obj at optimal: "<<objFn(std::get<0>(results), std::get<1>(results), std::get<2>(results), std::get<3>(results))<<std::endl;

    std::cout<<"obj at actual: "<<objFn(sig, lambda, muJ, sigJ)<<std::endl;
*/

}




