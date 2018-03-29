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
        return price/stock-transformStrike>0?transformStrike:0.0
    }
    /**
     * equation 2.9
     * */
    template<typename AssetValue, typename Strike, typename Discount>
    auto xJ(const AssetValue& stock, const Strike& strike, const Discount& discount){
        return log(strike*discount/stock);
    }
    constexpr int KNOT=0;
    constexpr int GAMMA=1;
    /**
     * Tuple is a std::tuple of x_j (knots) and O_j (gammas)
     * */
    template<typename T, typename Tuple>
    auto cspline(const num& x, const std::vector<Tuple>& knots_gamma){
        const auto minDiff=x-std::get<KNOT>(knots_gamma.front());//difference between given x and first knot
        const auto maxDiff=std::get<KNOT>(knots_gamma.back())-x;//difference between last knot and given x
        const auto span=std::get<KNOT>(knots_gamma.back())-std::get<KNOT>(knots_gamma.front());//span of knots
        const auto tripleMin=futilities::const_power(maxZeroOrNumber(minDiff), 3); //minDiff^3, if minDiff>0 else 0
        const auto tripleMax=futilities::const_power(maxZeroOrNumber(maxDiff), 3);//maxDiff^3, if minDiff>0 else 0
        int startFrom=1;
        int endFrom=1;
        return futilities::sum_subset(knots_gamma, startFrom, endFrom, [&](const auto& tuple, const auto& index){
            const auto lambda=(std::get<KNOT>(knots_gamma.back())-std::get<KNOT>(tuple))/span;
            const auto currDiff=x-std::get<KNOT>(tuple);//difference between given x and first knot
            const auto tripleCurr=futilities::const_power(maxZeroOrNumber(minDiff), 3); //minDiff^3, if minDiff>0 else 0
            return (tripleCurr+maxDiff*(lambda-1)-minDiff*lambda)*std::get<GAMMA>(knots_gamma[index+1]);
        })+std::get<GAMMA>(knots_gamma.front())+std::get<GAMMA>(knots_gamma[1])*x;
    }

    template<typename Strike, typename MarketPrice, typename AssetValue, typename Discount>
    auto generateOEstimate(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options, const AssetValue& stock, const Discount& discount){
        int numStrikes=strikes.size();
        const auto knots_gamma_tmp=futilities::for_each_parallel(0, numStrikes, [](const auto& index){
            return std::make_tuple(xJ(stock, strikes[index], discount), oJ(options[index], stock, strikes[index], discount));
        });
        return [knots_gamma = std::move(knots_gamma_tmp)](const auto& x){
            return cspline(x, knots_gamma);
        };
    }
    


}

#endif