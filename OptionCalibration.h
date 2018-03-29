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
    template<typename Strike, typename MarketPrice, typename AssetValue, typename Discount>
    auto generateFOEstimate(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options, const AssetValue& stock, const Discount& discount){
        int numStrikes=strikes.size();
        auto xV=futilities::for_each_parallel_copy(strikes, [&](const auto& strike, const auto& index){
            return xJ(stock, strike, discount); 
        });
        auto oV=futilities::for_each_parallel_copy(strikes, [&](const auto& strike, const auto& index){
            return oJ(options[index], stock, strike, discount);
        });


        auto knots_gamma_tmp=futilities::for_each_parallel(0, numStrikes, [&](const auto& index){
            const auto xj=xJ(stock, strikes[index], discount);
            return std::make_tuple(xJ(stock, strikes[index], discount), oJ(options[index], stock, strikes[index], discount));
        });
        const auto firstX=std::get<XJ>(knots_gamma_tmp.front());
        const auto lastX=std::get<XJ>(knots_gamma_tmp.back());
        const double incrementX=.05;//arbitrary
        knots_gamma_tmp.insert(
            knots_gamma_tmp.begin(), 
            std::make_tuple(firstX-incrementX,0.0)
        );
        knots_gamma_tmp.emplace_back(std::make_tuple(lastX+incrementX,0.0));

        return [knots_gamma = std::move(knots_gamma_tmp)](const auto& u){
            return fSpline(u, knots_gamma);
        };
    }
    


}

#endif