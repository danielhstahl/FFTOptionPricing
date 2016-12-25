#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include "RungeKutta.h"
#include "OptionPricing.h"
#include "fft.h"
#include <iostream>
#include "CharacteristicFunctions.h"

TEST_CASE("FSTS", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
    //auto initValue=50.0;
    
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
    auto myOptionsPrice=optionprice::FSTS(numX, xmax, 
        discount, 
        [&](const auto& logR){return payoff(K, logR, K);}, 
        BSCF
    );
    auto myXDomain=optionprice::getFSTSUnderlying(-xmax, xmax, K, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BS(myXDomain[i], discount(myXDomain[i]), K, sig*sqrt(T))).epsilon(.0001));
    }
}
TEST_CASE("FangOosterlee", "OptionPricing"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto K=50.0; 
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
    int numU=256;
    auto myXDomain=optionprice::getFangOostUnderlying(-xmax, xmax, K, numX);

    auto myOptionsPrice=optionprice::FangOost(numX, numU, xmax,K, discount, BSCF);
    
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){ 
        auto put=BS(myXDomain[i], discount, K, sig*sqrt(T))+K*discount-myXDomain[i];
        //std::cout<<myOptionsPrice[i+1]<<", "<<myXDomain[i+1]<<", "<<put<<", "<<put/myOptionsPrice[i+1]<<std::endl;
        REQUIRE(myOptionsPrice[i]==Approx(put).epsilon(.0001));
    }
}

TEST_CASE("CarrMadan", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=1.0;
    auto S0=5.0;
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
    auto myOptionsPrice=optionprice::CarrMadan(numX,  
        ada,
        S0, 
        discount, 
        [&](const auto& v, const auto& alpha){return optionprice::CallAug(v, alpha, BSCF);}
    );
    auto b=optionprice::getMaxK(ada);
    auto lambda=optionprice::getLambda(numX, b);
    auto myXDomain=optionprice::getCarrMadanStrikes(ada, S0, numX);
    int i=(int)numX*.3;
    int mxX=(int)numX*.7;
    for(;i<mxX; ++i){
        REQUIRE(myOptionsPrice[i]==Approx(BS(S0, discount, myXDomain[i], sig*sqrt(T))).epsilon(.0001));
    }
}



