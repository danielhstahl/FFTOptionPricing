#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include "RungeKutta.h"
#include "OptionPricing.h"
#include <chrono>
#include "fft.h"
#include <iostream>
#include "CharacteristicFunctions.h"

TEST_CASE("FSTSCall", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    
    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=[&](const auto& x){
        return exp(-r*T);
    };

    auto BS=[&](const auto &S0, const auto &discount, const auto &k, const auto &sigma){ //note that sigma includes sqrt(t) term so in vanilla BS sigma is equal to volatility*sqrt(T)
        if(sigma>0){
            double s=sqrt(2.0);
            auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
            return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma)/s)));
        }
        else{
            if(S0>k){
                return (S0-k)*discount;
            }
            else{
                return 0.0;
            }
        }
    };

    auto payoff=[&](const auto& initVal, const auto& logResult, const auto& K){
        auto assetValue=initVal*exp(logResult);
        return assetValue>K?assetValue-K:0.0;
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FSTS(numX, xmax, 
        discount, 
        [&](const auto& logR){return payoff(K, logR, K);}, 
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    auto myXDomain=optionprice::getFSTSUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BS(myXDomain[i], discount(myXDomain[i]), K, sig*sqrt(T))).epsilon(.0001));
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
    auto discount=[&](const auto& x){
        return exp(-r*T);
    };

    auto BS=[&](const auto &S0, const auto &discount, const auto &k, const auto &sigma){ //note that sigma includes sqrt(t) term so in vanilla BS sigma is equal to volatility*sqrt(T)
        if(sigma>0){
            double s=sqrt(2.0);
            auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
            return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma)/s)));
        }
        else{
            if(S0>k){
                return (S0-k)*discount;
            }
            else{
                return 0.0;
            }
        }
    };

    auto payoff=[&](const auto& initVal, const auto& logResult, const auto& K){
        auto assetValue=initVal*exp(logResult);
        return assetValue<K?K-assetValue:0.0;
    };
    int numX=pow(2, 10);
    double xmax=5.0;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FSTS(numX, xmax, 
        discount, 
        [&](const auto& logR){return payoff(K, logR, K);}, 
        BSCF
    );
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "FSTS time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    auto myXDomain=optionprice::getFSTSUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BS(myXDomain[i], discount(myXDomain[i]), K, sig*sqrt(T))-myXDomain[i]+K*discount(myXDomain[i])).epsilon(.0001));
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

    auto BS=[&](const auto &S0, const auto &discount, const auto &k, const auto &sigma){ //note that sigma includes sqrt(t) term so in vanilla BS sigma is equal to volatility*sqrt(T)
        if(sigma>0){
            double s=sqrt(2.0);
            auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
            return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma)/s)));
        }
        else{
            if(S0>k){
                return (S0-k)*discount;
            }
            else{
                return 0.0;
            }
        } 
    };
    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FangOostCall(S0, KArray, numU, discount, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticOption=BS(S0, discount, KArray[i], sig*sqrt(T));
        REQUIRE(myOptionsPrice[i]==Approx(analyticOption).epsilon(.0001));
    }
}
TEST_CASE("FangOosterleePut", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=50.0; 
    //auto initValue=50.0;
    

    auto BSCF=[&](const auto& u){
        return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
    };
    auto discount=exp(-r*T);

    auto BS=[&](const auto &S0, const auto &discount, const auto &k, const auto &sigma){ //note that sigma includes sqrt(t) term so in vanilla BS sigma is equal to volatility*sqrt(T)
        if(sigma>0){
            double s=sqrt(2.0);
            auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
            return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma)/s)));
        }
        else{
            if(S0>k){
                return (S0-k)*discount;
            }
            else{
                return 0.0;
            }
        } 
    };
    double xmax=5.0;
    int numX=pow(2, 10);
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto KArray=optionprice::getFangOostStrike(-xmax, xmax, S0, numX);
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FangOostPut(S0, KArray, numU, discount, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;

    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto analyticOption=BS(S0, discount, KArray[i], sig*sqrt(T))+KArray[i]*discount-S0;
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
    auto BS=[&](const auto &S0, const auto &discount, const auto &k, const auto &sigma){ //note that sigma includes sqrt(t) term so in vanilla BS sigma is equal to volatility*sqrt(T)
        if(sigma>0){
            double s=sqrt(2.0);
            auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
            return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma)/s)));
        }
        else{
            if(S0>k){
                return (S0-k)*discount;
            }
            else{
                return 0.0;
            }
        }
    };
    int numX=pow(2, 10);
    auto ada=.25;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::CarrMadanCall(numX,  
        ada,
        S0, 
        discount, 
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
        REQUIRE(myOptionsPrice[i]==Approx(BS(S0, discount, myXDomain[i], sig*sqrt(T))).epsilon(.0001));
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
    auto BS=[&](const auto &S0, const auto &discount, const auto &k, const auto &sigma){ //note that sigma includes sqrt(t) term so in vanilla BS sigma is equal to volatility*sqrt(T)
        if(sigma>0){
            double s=sqrt(2.0);
            auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
            return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma)/s)));
        }
        else{
            if(S0>k){
                return (S0-k)*discount;
            }
            else{
                return 0.0;
            }
        }
    };
    int numX=pow(2, 10);
    auto ada=.25;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        discount, 
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
        REQUIRE(myOptionsPrice[i]==Approx(BS(S0, discount, myXDomain[i], sig*sqrt(T))-S0+myXDomain[i]*discount).epsilon(.0001));
    }
}



