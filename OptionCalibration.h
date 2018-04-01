#ifndef OPTIONCALIBRATION_H
#define OPTIONCALIBRATION_H
#include <complex>
#include <vector>
#include <algorithm>
#include "FunctionalUtilities.h"
#include "fft.h"
#include <tuple>
#include "Newton.h"
/**
 * Optimization specific
 * to Levy processes and
 * the inverse problem for 
 * option prices.
 * Based off the following 
 * work:
 * Option calibration of 
 * exponential Levy models
 * by Belomestny and ReiB
 * */
namespace optioncal{
    /**
     * equation 2.8
     * */
    template<typename MarketPrice, typename AssetValue, typename Strike, typename Discount>
    auto oJ(const MarketPrice& price, const AssetValue& stock, const Strike& strike,  const Discount& discount){
        auto transformStrike=1.0-strike*discount/stock;
        auto transformStrikePlus=transformStrike>0?transformStrike:0.0;
        return price/stock-transformStrikePlus;
    }
    /**
     * equation 2.9
     * */
    template<typename AssetValue, typename Strike, typename Discount>
    auto xJ(const AssetValue& stock, const Strike& strike, const Discount& discount){
        return log(strike*discount/stock);
    }
    template<typename T>
    auto maxZeroOrNumber(const T& num){
        return num>0.0?num:T(0.0);
    }
    constexpr int XJ=0;
    constexpr int OJ=1;



    /**
    u is real (not complex)
    */
    template<typename U, typename Tuple>
    auto fSpline(const U& u, const std::vector<Tuple>& knots_gamma){
        const int startFrom=1;
        const int endFrom=1;
        return futilities::sum_subset(knots_gamma, startFrom, endFrom, [&](const auto& tuple, const auto& index){
            const auto currX=std::get<XJ>(tuple);
            const auto prevX=std::get<XJ>(knots_gamma[index-1]);
            const auto nextX=std::get<XJ>(knots_gamma[index+1]);
            const auto currO=std::get<OJ>(tuple);

            const auto currExp=exp(std::complex<U>(0.0, currX*u));
            const auto prevExp=exp(std::complex<U>(0.0, prevX*u));
            const auto nextExp=exp(std::complex<U>(0.0, nextX*u));

            //equation 3.10
            const auto retVal=currO*((currExp-prevExp)/(currX-prevX)-(nextExp-currExp)/(nextX-currX))/(u*u);
            //find the "j0" and "j0-1"
            return currX<0&&nextX>=0?retVal+(1.0+(nextExp*currX-currExp*nextX)/(nextX-currX))/(u*u):retVal;
        });
    }

    /*auto bound(double maxStrike){
        return log(maxStrike);
    }*/
    
    /**STRIKES NEEDS TO BE IN ORDER!!  The strikes and options are passed as values since we need to modify them*/
   /* auto generateFOEstimate(std::vector<double> strikes, std::vector<double> options, const double& stock, const double& maxStrike){
        int numStrikes=strikes.size();
        strikes=futilities::for_each_parallel(std::move(strikes), [&](const auto& strike, const auto& index){
            return xJ(stock, strike, discount); 
        });
        options=futilities::for_each_parallel(std::move(options), [&](const auto& option, const auto& index){
            return oJ(option, stock, strikes[index], discount);
        });


        options.insert(options.begin(), stock); //option price is stock price at k=0
        options.push_back(0.0);//option price is zero for K->infty
        
        strikes=futilities::for_each_parallel(std::move(strikes), [](const auto& k, const auto& index){
            return log(k);
        });
        double lStrike=bound(maxStrike);
        strikes.insert(strikes.begin(), -lStrike);
        strikes.push_back(lStrike);

        //tk::spline sTmp;
        //sTmp.set_points(strikes, options); 
        return [s = std::move(sTmp)](const auto& callPrice){
            return s(callPrice);
        };
    }*/

    template<typename Strike, typename MarketPrice, typename AssetValue, typename Discount>
    auto generateFOEstimate(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options, const AssetValue& stock, const Discount& discount, const Strike& maxStrike){
        int numStrikes=strikes.size();

        auto knots_gamma_tmp=futilities::for_each_parallel(0, numStrikes, [&](const auto& index){
            const auto xj=xJ(stock, strikes[index], discount);
            return std::make_tuple(xJ(stock, strikes[index], discount), oJ(options[index], stock, strikes[index], discount));
        });

        knots_gamma_tmp.insert(
            knots_gamma_tmp.begin(), 
            std::make_tuple(
                xJ(stock, 1/maxStrike, discount),
                0.0
            )
        );
        knots_gamma_tmp.emplace_back(std::make_tuple(xJ(stock, maxStrike, discount),0.0));

        return [knots_gamma = std::move(knots_gamma_tmp)](const auto& u){
            const auto vi=std::complex<double>(0, u);
            //equation 3.1
            return log(1.0+vi*(1.0+vi)*fSpline(u, knots_gamma));
        };
    }




    /*constexpr int DK=0;
    constexpr int DU=1;
    auto getDKandDU(int N, double maxStrike){
        double lStrike=bound(maxStrike);
        return std::make_tuple(2*lStrike/N, M_PI/lStrike);
    }

    template<typename IFS, typename Discount>
    auto cfHat(
        const IFS& instSpline, 
        const Discount& discount, 
        double maxStrike,
        double dk, double du, int N
    ){
        double lStrike=bound(maxStrike);
        return futilities::for_each_parallel(
            ifft(futilities::for_each_parallel(0, N, [&](const auto& index){
                auto pm=index%2==0?-1.0:1.0; //simpson's rule
                auto currK=-lStrike+dk*index;
                return std::complex<double>(instSpline(currK)*exp(-dk*index)*dk*(3.0+pm)/3.0, 0.0);
            })),
            [&](const auto& xn, const auto& index){
                auto currU=du*index;
                auto v=std::complex<double>(-1.0, currU);
                return -v*v*xn*discount*exp(lStrike*(1.0-std::complex<double>(0.0, currU)));
            }
        );
    }*/


    template<typename PhiHatFn, typename LogCfFN, typename DiscreteU>
    auto getObjFn(PhiHatFn&& phiHatFntmp, LogCfFN&& cfFntmp, DiscreteU&& uArraytmp){
        return [phiHatFn=std::move(phiHatFntmp), cfFn=std::move(cfFntmp), uArray=std::move(uArraytmp)](const auto&... params){
            return futilities::sum(uArray, [&](const auto& u, const auto& index){
                //u-i comes from equation 3.1
                return std::norm(
                    phiHatFn(u)-cfFn(u-std::complex<double>(0, 1.0), params...)
                );
            })/uArray.size();
        };
    }
    /**fn is the result from getObjFn*/
    template<typename FN, typename ...Params>
    auto calibrate(const FN& fn, const Params&... params){
        const int maxIter=50;
        const double prec=.00001;
        const double peterb=.000001;
        const double alpha=.0001; //needs a very small step or it goes off to no where
        return newton::gradientDescentApprox(fn, maxIter, prec, peterb, alpha, params...);
    }

        


}

#endif