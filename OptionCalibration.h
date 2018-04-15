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
//#include "optim.hpp"


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
    template<typename AssetValue, typename Strike, typename Discount>
    auto kX(const AssetValue& stock, const Strike& x, const Discount& discount){
        return stock*exp(x)/discount;
    }
    template<typename T>
    auto maxZeroOrNumber(const T& num){
        return num>0.0?num:T(0.0);
    }
    constexpr int XJ=0;
    constexpr int OJ=1;
    constexpr auto cmp=std::complex<double>(0, 1);


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

    /*double BSCall(double S0, double discount, double k, double sigma, double T){
        double s=sqrt(2.0);
        double sigmasqrt=sqrt(T)*sigma;
        auto d1=log(S0/(discount*k))/(sigmasqrt)+sigmasqrt*.5;
        return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigmasqrt)/s)));
    }*/

    template<typename Strike, typename MarketPrice, typename AssetValue, typename Discount>
    auto generateFOEstimateSp(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options, const AssetValue& stock, const Discount& discount, const Strike& maxStrike){
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

    auto getDU(int N, double uMax){
        return 2.0*uMax/N;
    }
    auto getDX(int N, double xMin,double xMax){
        return (xMax-xMin)/(N-1);
    }
    auto getUMax(int N, double xMin, double xMax){
        return M_PI*(N-1)/(xMax-xMin);
    }
    
    template<typename PhiHatFn, typename LogCfFN, typename DiscreteU, typename Strike, typename AssetValue, typename Discount>
    auto getObjFnSpline(PhiHatFn&& phiHatFntmp, LogCfFN&& cfFntmp, int N, const Strike& minStrike, const Strike& maxStrike, const AssetValue& stock , const Discount& discount){
        auto xMin=xJ(stock, minStrike, discount);
        auto xMax=xJ(stock, maxStrike, discount);
        return [phiHatFn=std::move(phiHatFntmp), cfFn=std::move(cfFntmp), N, xM=std::move(xMin), xP=std::move(xMax)](const auto&... params){
            const auto uMax=getUMax(N, xM, xP);
            const auto du=getDU(N, uMax);
            auto uArray=futilities::for_each_parallel(0, N, [&](const auto& index){
                return index*du-uMax;
            }); 
            //return 
            //taking average of u over domain and then returning the norm of the difference
            return std::norm(
                (
                    futilities::sum(phiHatFn(N), [&](const auto& u, const auto& index){
                        return u; 
                    })-
                    futilities::sum(uArray, [&](const auto& u, const auto& index){
                        return cfFn(u*cmp, params...);
                    })
                )/(double)uArray.size()           
            );
        };
    }

    template<typename ArrayU, typename Fn, typename ApplyEachSum, typename T>
    auto DFT(const std::vector<ArrayU>& uArray, Fn&& fn, ApplyEachSum&& fnU, const T& xMin, const T& xMax, int N){
        T dx=(xMax-xMin)/(N-1);
        int uSize=uArray.size();
        return futilities::for_each_parallel(0, uSize, [&](const auto& uIndex){
            auto u=uArray[uIndex];
            return fnU(u, futilities::sum(0, N, [&](const auto& index){
                auto simpson=(index==0||index==(N-1))?1.0:(index%2==0?2.0:4.0);
                auto x=xMin+dx*index;
                return exp(cmp*u*x)*fn(u, x, index)*simpson*dx/3.0;
            }));
        });
        
    }

    template<typename StrikeArray, typename Strike, typename MarketPrice, typename AssetPrice, typename Discount>
    auto generateFOEstimate(const StrikeArray& strikes, const std::vector<MarketPrice>& options, const Discount& discount, const AssetPrice& stock, const Strike& minStrike, const Strike& maxStrike){
        const int numStrikes=strikes.size();
        StrikeArray paddedStrikes(numStrikes+2);
        paddedStrikes=futilities::for_each_parallel_subset(std::move(paddedStrikes), 1, 1, [&](const auto& v, const auto& i){
            return strikes[i-1]/stock;
        });
        paddedStrikes[0]=minStrike/stock;
        paddedStrikes[numStrikes+1]=maxStrike/stock;
        std::vector<MarketPrice> optionPrices(numStrikes+2);
        optionPrices=futilities::for_each_parallel_subset(std::move(optionPrices), 1, 1, [&](const auto& v, const auto& i){
            return options[i-1]/stock;
        });
        optionPrices[0]=(stock-minStrike*discount)/stock; //close to just being the stock price
        optionPrices[numStrikes+1]=0.0; //just 0 at the end

        tk::spline s;
        s.set_points(paddedStrikes, optionPrices);
        /*for(int i=0; i<128; ++i){
            double strike=-5.0+i*(8.0)/127.0;
            std::cout<<"strike: "<<exp(strike)<<", option: "<<s(exp(strike))<<std::endl;
        }*/
        return [
            spline=std::move(s), 
            minStrike=std::move(paddedStrikes.front()), 
            maxStrike=std::move(paddedStrikes.back()),
            discount=std::move(discount)
        ](int N, const auto& uArray){
            const auto xMin=log(discount*minStrike);
            const auto xMax=log(discount*maxStrike);
            auto valOrZero=[](const auto& v){
                return v>0.0?v:0.0;
            };
            return DFT(uArray, [&](const auto& u, const auto& x, const auto& index){
                const auto expX=exp(x);
                const auto strike=expX/discount;
                const auto offset=1.0-expX;
                const auto optionPrice=spline(strike);
                //const auto optionPrice=BSCall(1.0, discount, strike, .3, 1.0);
                return valOrZero(optionPrice)-valOrZero(offset);
            }, [&](const auto& u, const auto& value){
                const auto front=u*cmp*(1.0+u*cmp);
                return log(1.0+front*value);
            }, xMin, xMax, N);
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


    template<typename PhiHat, typename LogCfFN, typename DiscreteU>
    auto getObjFn(PhiHat&& phiHattmp, LogCfFN&& cfFntmp, DiscreteU&& uArraytmp){
        return [phiHat=std::move(phiHattmp), cfFn=std::move(cfFntmp), uArray=std::move(uArraytmp)](const auto&... params){
            return futilities::sum(uArray, [&](const auto& u, const auto& index){
                return std::norm(
                    phiHat[index]-cfFn(std::complex<double>(1.0, u), params...)
                )/(double)uArray.size();         
            });
        };
    }
    template<typename PhiHat, typename LogCfFN, typename DiscreteU>
    auto getObjFn_arr(PhiHat&& phiHattmp, LogCfFN&& cfFntmp, DiscreteU&& uArraytmp){
        return [phiHat=std::move(phiHattmp), cfFn=std::move(cfFntmp), uArray=std::move(uArraytmp)](const auto& arrayParam){
            return futilities::sum(uArray, [&](const auto& u, const auto& index){
                return std::norm(
                    phiHat[index]-cfFn(std::complex<double>(1.0, u), arrayParam)
                )/(double)uArray.size();         
            });
        };
    }
    /**fn is the result from getObjFn*/
    template<typename FN, typename ...Params>
    auto calibrate(const FN& fn, const Params&... params){
        const int maxIter=500;
        const double prec=.00001; 
        const double peterb=.000001;
        static const std::size_t value = sizeof...(Params);
        const double alpha=.01; //needs a very small step or it goes off to no where
        return newton::gradientDescentApprox(fn, maxIter, prec, peterb, alpha, params...);
    }



}

#endif