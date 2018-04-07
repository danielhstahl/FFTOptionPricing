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
    knots_gamma MUST be in order (least to greatest in x)
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

            const auto currExp=exp(std::complex<double>(0.0, currX)*u);
            const auto prevExp=exp(std::complex<double>(0.0, prevX)*u);
            const auto nextExp=exp(std::complex<double>(0.0, nextX)*u);

            //equation 3.10
            const auto retVal=currO*((currExp-prevExp)/(currX-prevX)-(nextExp-currExp)/(nextX-currX))/(u*u);
            //find the "j0" and "j0-1"
            return currX<0&&nextX>=0?retVal+(1.0+(nextExp*currX-currExp*nextX)/(nextX-currX))/(u*u):retVal;
        });
    }

 
    template<typename Strike, typename MarketPrice, typename AssetValue, typename Discount>
    auto generateFOEstimate(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options, const AssetValue& stock, const Discount& discount, const Strike& maxStrike){
        int numStrikes=strikes.size();
        auto knots_gamma_tmp=std::vector<std::tuple<Strike, AssetValue> >(numStrikes+2);
        int lengthFromEdge=1;
        knots_gamma_tmp=futilities::for_each_parallel_subset(std::move(knots_gamma_tmp), lengthFromEdge, lengthFromEdge, [&](const auto& v, const auto& index){
            const auto pIndex=index-lengthFromEdge;
            const auto xj=xJ(stock, strikes[pIndex], discount);
            return std::make_tuple(xJ(stock, strikes[pIndex], discount), oJ(options[pIndex], stock, strikes[pIndex], discount));
        });

        knots_gamma_tmp[0]=std::make_tuple(
            xJ(stock, stock/maxStrike, discount),
            0.0
        );
        knots_gamma_tmp[numStrikes+lengthFromEdge]=std::make_tuple(
            xJ(stock, stock*maxStrike, discount),
            0.0
        );

        return [knots_gamma = std::move(knots_gamma_tmp)](const auto& u){
            const auto vi=std::complex<double>(0.0, 1.0)*u;
            //equation 3.1
            return u==0?vi:log(1.0+vi*(1.0+vi)*fSpline(u, knots_gamma));
        };
    }

    template<typename PhiHatFn, typename LogCfFN, typename DiscreteU>
    auto getObjFn(PhiHatFn&& phiHatFntmp, LogCfFN&& cfFntmp, DiscreteU&& uArraytmp){
        return [phiHatFn=std::move(phiHatFntmp), cfFn=std::move(cfFntmp), uArray=std::move(uArraytmp)](const auto&... params){
            return futilities::sum(uArray, [&](const auto& u, const auto& index){
                //u-i comes from equation 3.1
                return std::norm(
                    phiHatFn(u)-cfFn(u*std::complex<double>(0, 1.0)+1.0, params...)
                );
            })/uArray.size();
        };
    }
    /**fn is the result from getObjFn*/
    template<typename FN, typename ...Params>
    auto calibrate(const FN& fn, const Params&... params){
        const int maxIter=500;
        const double prec=.00001; 
        const double peterb=.000001;
        const double alpha=.01; //needs a very small step or it goes off to no where
        return newton::gradientDescentApprox(fn, maxIter, prec, peterb, alpha, params...);
    }

}

#endif