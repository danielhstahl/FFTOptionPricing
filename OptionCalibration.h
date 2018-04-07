#ifndef OPTIONCALIBRATION_H
#define OPTIONCALIBRATION_H
#include <complex>
#include <vector>
#include <algorithm>
#include "FunctionalUtilities.h"
#include "fft.h"
#include <tuple>
#include "spline.h"
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
 * and the updates from Sohl
 * and Trabs
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
            //const auto currExp=exp(currX*u);
            const auto prevExp=exp(std::complex<double>(0.0, prevX)*u);
            //const auto prevExp=exp(prevX*u);
            const auto nextExp=exp(std::complex<double>(0.0, nextX)*u);
            //const auto nextExp=exp(nextX*u);

            //equation 3.10 
            const auto retVal=currO*((currExp-prevExp)/(currX-prevX)-(nextExp-currExp)/(nextX-currX))/(u*u);
            //find the "j0-1" and "j0"
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
            const auto uPlusi=std::complex<double>(0.0, 1.0)+u;
            //equation 3.1...note that we are solving for the single argument "u" into the analytical CF
            return u==0?0.0:log(1.0-uPlusi*u*fSpline(uPlusi, knots_gamma));
        };
    }

    auto getUMax(int N, double xMax){
        return M_PI*N/(2*xMax);
    }
    auto getDU(int N, double uMax){
        return 2.0*uMax/N;
    }


    template<typename Strike, typename MarketPrice, typename AssetValue, typename Discount>
    auto generateFOEstimateSpline(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options, const AssetValue& stock, const Discount& discount, const Strike& maxStrike){
        int numStrikes=strikes.size();
        int lengthFromEdge=1;
        std::vector<MarketPrice> O(numStrikes+2, 0.0);
        std::vector<Strike> X(numStrikes+2);
        O=futilities::for_each_parallel_subset(std::move(O), lengthFromEdge, lengthFromEdge, [&](const auto& v, const auto& index){
            const auto pIndex=index-lengthFromEdge;
            return oJ(options[pIndex], stock, strikes[pIndex], discount);
        });
        X=futilities::for_each_parallel_subset(std::move(X), lengthFromEdge, lengthFromEdge, [&](const auto& v, const auto& index){
            const auto pIndex=index-lengthFromEdge;
            return xJ(stock, strikes[pIndex], discount);
        });
        X[0]=xJ(stock, stock/maxStrike, discount);
        X[numStrikes+lengthFromEdge]=xJ(stock, maxStrike, discount);

        tk::spline s;
        s.set_points(X,O); 

        return [spline=std::move(s)](const auto& N, const auto& xMax){
            const auto uMax=getUMax(N, xMax);
            const auto dx=2.0*xMax/N;
            const auto du=getDU(N, uMax);
            //dx=2a/N, du=pi/a=2pi/(dx*N)
            return futilities::for_each_parallel(
                ifft(futilities::for_each_parallel(0, N, [&](const auto& index){
                    auto pm=index%2==0?-1.0:1.0; //simpson's rule
                    auto currX=dx*index-xMax;
                    auto xOffset=exp(-(std::complex<double>(0.0, uMax)+1.0)*currX);
                    auto abOffset=exp(xMax*uMax*std::complex<double>(0, 1));
                    return abOffset*xOffset*spline(currX)*dx*(3.0+pm)/3.0;
                })),
                [&](const auto& xn, const auto& index){
                    auto u=du*index-uMax;
                    const auto splineResult=exp(-std::complex<double>(0, 1)*u*xMax)*xn;
                    const auto uPlusi=std::complex<double>(0.0, 1.0)+u;
                    //equation 3.1...note that we are solving for the single argument "u" into the analytical CF
                    return u==0?0.0:log(1.0-uPlusi*u*splineResult);
                }
            );
        };
    }
    template<typename PhiHatFn, typename LogCfFN, typename DiscreteU>
    auto getObjFnSpline(PhiHatFn&& phiHatFntmp, LogCfFN&& cfFntmp, int N, double xMax){
        return [phiHatFn=std::move(phiHatFntmp), cfFn=std::move(cfFntmp), N, xMax](const auto&... params){
            double uMax=getUMax(N, xMax);
            const auto du=getDU(N, uMax);
            auto uArray=futilities::for_each_parallel(0, N, [&](const auto& index){
                return index*du-uMax;
            }); 
            //return 
            //taking average of u over domain and then returning the norm of the difference
            return std::norm(
                (
                    futilities::sum(phiHatFn(N, xMax), [&](const auto& u, const auto& index){
                        return u; 
                    })-
                    futilities::sum(uArray, [&](const auto& u, const auto& index){
                        return cfFn(u*std::complex<double>(0, 1.0), params...);
                    })
                )/(double)uArray.size()           
            );
        };
    }

/**

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
}*/


    template<typename PhiHatFn, typename LogCfFN, typename DiscreteU>
    auto getObjFn(PhiHatFn&& phiHatFntmp, LogCfFN&& cfFntmp, DiscreteU&& uArraytmp){
        return [phiHatFn=std::move(phiHatFntmp), cfFn=std::move(cfFntmp), uArray=std::move(uArraytmp)](const auto&... params){
            
            //return 
            //taking average of u over domain and then returning the norm of the difference
            return std::norm(
                (
                    futilities::sum(uArray, [&](const auto& u, const auto& index){
                        return phiHatFn(u); 
                    })-
                    futilities::sum(uArray, [&](const auto& u, const auto& index){
                        return cfFn(u*std::complex<double>(0, 1.0), params...);
                    })
                )/(double)uArray.size()           
            );
        };
    }
    /**fn is the result from getObjFn*/
    template<typename FN, typename ...Params>
    auto calibrate(const FN& fn, const Params&... params){
        const int maxIter=500;
        const double prec=.00001; 
        const double peterb=.000001;
        static const std::size_t value = sizeof...(Params);
        const double alpha=.01;//*value; //needs a very small step or it goes off to no where
        return newton::gradientDescentApprox(fn, maxIter, prec, peterb, alpha, params...);
    }

}

#endif