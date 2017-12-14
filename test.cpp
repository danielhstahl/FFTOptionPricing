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



TEST_CASE("FangOosterleeCallFewK", "[OptionPricing]"){
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
    int numU=64; //This should run extremely quickly since only have 5 strikes
    auto KArray=std::vector<double>(5, 0);
    KArray[4]=.3;
    KArray[3]=40;//strike 40
    KArray[2]=50;//strike 50
    KArray[1]=60;//strike 60
    KArray[0]=7500;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FangOostCall(S0, KArray, numU, discount, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    for(int i=1;i<4; ++i){ 
        auto analyticOption=BS(S0, discount, KArray[i], sig*sqrt(T));
        REQUIRE(myOptionsPrice[i]==Approx(analyticOption).epsilon(.0001));
    }
}
TEST_CASE("FangOosterleeCallFewKLowT", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=.25;
    auto S0=50.0; 
   

    auto BSCF=[&](const auto& u){
        //return chfunctions::gaussCF(u, (r-sig*sig*.5)*T, sig*sqrt(T));
        return exp(chfunctions::gaussLogCF(u, r-sig*sig*.5, sig)*T);
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
    int numU=256; //This should run extremely quickly since only have 5 strikes, the higher U is to compensate for low t
    auto KArray=std::vector<double>(5, 0);
    KArray[4]=.3;
    KArray[3]=40;//strike 40
    KArray[2]=50;//strike 50
    KArray[1]=60;//strike 60
    KArray[0]=7500;
    auto started = std::chrono::high_resolution_clock::now();
    auto myOptionsPrice=optionprice::FangOostCall(S0, KArray, numU, discount, BSCF);
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Fang Oosterlee time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    for(int i=1;i<4; ++i){ 
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
TEST_CASE("CarrMadanCallLowT", "[OptionPricing]"){
    auto r=.05;
    auto sig=.3;
    auto T=.25;
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
    auto discount=exp(-r*T);

    auto BSCF=[&](const auto& u){
        return exp(chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T);
    };
    
    int numU=64; //this will be slightly slower than Carr madan or FSTS....since carr and madan run in nlogn where n=1024 and 64 is roughly 10 times larger than log(1024)
    auto myOptionsPrice=optionprice::FangOostCall(S0, KArray, numU, discount, BSCF);
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
    auto discount=exp(-r*T);

    auto BSCF=[&](const auto& u){
        return exp(chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T);
    };
    
    int numU=64; //
    auto myOptionsPrice=optionprice::FangOostCall(S0, KArray, numU, discount, BSCF);
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
    
    int numU=256; //higher due to low T
    auto myOptionsPrice=optionprice::FangOostCall(S0, KArray, numU, discount, BSCF);
    auto myReference=BS(S0, discount, KArray[1], sig*sqrt(T));
    REQUIRE(myOptionsPrice[1]==Approx(myReference).epsilon(.0001));


}


TEST_CASE("FangOosterleeCallCGMYWithVol", "[OptionPricing]"){
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
    auto myOptionsPrice=optionprice::FangOostPut(S0, KArray, numU, discount, BSCF);
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
    auto discount=exp(-r*T);

    auto BSCF=[&](const auto& u){
        return exp(chfunctions::cgmyLogRNCF(u, C, G, M, Y, r, sig)*T);
    };
    
    
    int numX=pow(2, 10); //this is high...this seems to be a computationally tricky problem
    auto ada=.25;
    auto myOptionsPrice=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        discount, 
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
    auto discount=exp(-r*T);
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
        discount, 
        CFCorr
    );
    auto myOptionsPriceStoch=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        discount, 
        CFStoch
    );
    auto myOptionsPriceBase=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        discount, 
        CFBase
    );
    auto myOptionsPriceCorrBM=optionprice::CarrMadanPut(numX,  
        ada,
        S0, 
        discount, 
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

