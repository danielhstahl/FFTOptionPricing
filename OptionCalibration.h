#ifndef OPTIONCALIBRATION_H
#define OPTIONCALIBRATION_H
#include <complex>
#include <vector>
#include <algorithm>
#include "FunctionalUtilities.h"
#include "fft.h"
#include <tuple>
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
        return price/stock-transformStrike>0?transformStrike:0.0;
    }
    /**
     * equation 2.9
     * */
    template<typename AssetValue, typename Strike, typename Discount>
    auto xJ(const AssetValue& stock, const Strike& strike, const Discount& discount){
        return log(strike*discount/stock);
    }
    constexpr int XJ=0;
    constexpr int OJ=1;
    /**
     * Tuple is a std::tuple of x_j (knots) and O_j (gammas)
     * U is not complex
     * */
    template<typename U, typename Tuple>
    auto fSpline(const U& u, const std::vector<Tuple>& knots_gamma){
        
        const int startFrom=1;
        const int endFrom=1;


        return futilities::sum_subset(knots_gamma, startFrom, endFrom, [&](const auto& tuple, const auto& index){
            const auto currX=std::get<XJ>(tuple);
            const auto prevX=std::get<XJ>(knots_gamma[index-1]);
            const auto nextX=std::get<XJ>(knots_gamma[index+1]);
            const auto currO=std::get<OJ>(tuple);

            const auto currExp=exp(std::complex<U>(0.0, currX));
            const auto prevExp=exp(std::complex<U>(0.0, prevX));
            const auto nextExp=exp(std::complex<U>(0.0, nextX));

            //equation 3.10
            const auto retVal=currO*((currExp-prevExp)/(currX-prevX)-(nextExp-currExp)/(nextX-currX))/(u*u);
            //find the "j0" and "j0-1"
            return currX<0&&nextX>=0?retVal+(1.0+(nextExp*currX-currExp*nextX)/(nextX-currX))/(u*u):retVal;
        });
    }

   /* template<typename Strike, typename MarketPrice, typename AssetValue, typename Discount>
    auto generateFOEstimate(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options, const AssetValue& stock, const Discount& discount){
        int numStrikes=strikes.size();
        const auto knots_gamma_tmp=futilities::for_each_parallel(0, numStrikes, [](const auto& index){
            return std::make_tuple(xJ(stock, strikes[index], discount), oJ(options[index], stock, strikes[index], discount));
        });
        return [knots_gamma = std::move(knots_gamma_tmp)](const auto& x){
            return cspline(x, knots_gamma);
        };
    }*/
    
    /**STRIKES NEEDS TO BE IN ORDER!!*/
    template<typename Strike, typename MarketPrice>
    auto generateFOEstimate(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options, const Asset& stock){
        int numStrikes=strikes.size();
        auto knots_gamma_tmp=futilities::for_each_parallel(0, numStrikes, [&](const auto& index){
            return std::make_tuple(strikes[index], options[index]);
        });
        knots_gamma_tmp.insert(
            knots_gamma_tmp.begin(), 
            std::make_tuple(0.0,stock)//at k=0, the call equals the stock
        );
        constexpr double upperBound=100;//arbitrary
        knots_gamma_tmp.emplace_back(std::make_tuple(options.back()*upperBound,0.0)); //at large values of strike, option is zero

        return [knots_gamma = std::move(knots_gamma_tmp)](const auto& callPrice){
            return fSpline(callPrice, knots_gamma);
        };
    }

    template<typename IFS, typename Discount>
    auto getCFFromMarketData(const IFS& instSpline, const std::vector<MarketPrice>& options, const Discount& discount, int N){
        constexpr double upperBound=100;//arbitrary...this may be too high
        const maxStrike=options.back()*upperBound;
        const double dk=maxStrike/N;
        const double du=2.0*M_PI/maxStrike; ///may want to compute this outside thsi function since this will be needed to discretize the CF
        ///this should be min and max of call options.  The result iof the fft should be an approximate CF
        return futilities::for_each_parallel(
            ifft(futilities::for_each_parallel(0, N, [&](const auto& index){
                auto pm=index%2==0?-1.0:1.0; //simpson's rule
                auto currK=dk*index;
                return instSpline(currK)*currK*currK*(3.0+pm)/3.0;
            })),
            [&](const auto& xn, const auto& index){
                auto currU=du*index;
                return currU*currU*xn;
            }
        );
    }
    


}

#endif