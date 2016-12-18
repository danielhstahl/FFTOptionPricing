#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include "RungeKutta.h"
#include "OptionPricing.h"
#include "CharacteristicFunctions.h"

TEST_CASE("BlackScholes", "[Functional]"){
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
    int numX=pow(2, 8);
    auto myOptionsPrice=optionprice::FSTS(numX, -5.0, 5.0, discount, [&](const auto& logR){return payoff(K, logR, K);}, BSCF);
    auto myXDomain=futilities::for_emplace_back(-5.0, 5.0, numX, [&](const auto& val){
        return K*exp(val);
    });
    futilities::for_each_parallel(myXDomain, [&](const auto& val, const auto& index){
        REQUIRE(myOptionsPrice[index]==BS(val, discount(val), K, sig));
        return 0;
    });

    
}

TEST_CASE("BlackScholes", "[Functional]"){
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
    int numX=pow(2, 8);
    auto myOptionsPrice=optionprice::FSTS(numX, -5.0, 5.0, discount, [&](const auto& logR){return payoff(K, logR, K);}, BSCF);
    auto myXDomain=futilities::for_emplace_back(-5.0, 5.0, numX, [&](const auto& val){
        return K*exp(val);
    });
    futilities::for_each_parallel(myXDomain, [&](const auto& val, const auto& index){
        REQUIRE(myOptionsPrice[index]==BS(val, discount(val), K, sig));
        return 0;
    });

    
}