#ifndef OPTIONCALIBRATION_H
#define OPTIONCALIBRATION_H
#include <complex>
#include <vector>
#include <algorithm>
#include "FunctionalUtilities.h"
#include "fft.h"
#include <tuple>
#include "monotone_spline.h"
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
            const auto prevExp=exp(std::complex<double>(0.0, prevX)*u);
            const auto nextExp=exp(std::complex<double>(0.0, nextX)*u);

            //equation 3.10 
            const auto retVal=currO*((currExp-prevExp)/(currX-prevX)-(nextExp-currExp)/(nextX-currX))/(u*u);
            //find the "j0-1" and "j0"
            return currX<0&&nextX>=0?retVal+(1.0+(nextExp*currX-currExp*nextX)/(nextX-currX))/(u*u):retVal;
        });
    }

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

    template<typename T>
    auto transformPrice(const T& p, const T& v){
        return p/v;
    }

    template<typename Stock, typename V, typename FN>
    auto transformPrices(const std::vector<V>& arr, const Stock& asset, const V& minV, const V& maxV, FN&& transform){
        const int arrLength=arr.size();
        std::vector<V> paddedArr(arrLength+2);
        paddedArr=futilities::for_each_parallel_subset(std::move(paddedArr), 1, 1, [&](const auto& v, const auto& index){
            return transform(arr[index-1], asset, index);
        });
        paddedArr.front()=transform(minV, asset, 0);
        paddedArr.back()=transform(maxV, asset, arrLength);
        return std::move(paddedArr);
    }
    template<typename Stock, typename V>
    auto transformPrices(const std::vector<V>& arr, const Stock& asset, const V& minV, const V& maxV){
        return transformPrices(arr, asset, minV, maxV, [](const auto& p, const auto& v, const auto& index){
           return transformPrice(p, v);
        });
    }

    template<typename Arr, typename FN1, typename FN2, typename FN3>
    auto filter(const Arr& arr, FN1&& cmp, FN2&& fv1, FN3&& fv2){
        Arr v1;
        Arr v2;
        for(int i=0; i<arr.size();++i){
            const auto x=arr[i];
            if(cmp(x, i)){
                v1.emplace_back(fv1(x, i));
            }else{
                v2.emplace_back(fv2(x, i));
            }
        }
        return std::make_tuple(v1, v2);
    }

    template< typename Strike, typename MarketPrice, typename AssetPrice, typename Discount>
    auto getOptionSpline(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options,  const AssetPrice& stock, const Discount& discount, const Strike& minStrike, const Strike& maxStrike){
        auto paddedStrikes=transformPrices(strikes, stock, minStrike, maxStrike);
        const auto minOption=stock-minStrike*discount;
        const auto maxOption=0.0000001;
        auto valOrZero=[](const auto& v){
            return v>0.0?v:0.0;
        };
        auto optionPrices=transformPrices(options, stock, minOption, maxOption);
       
        const auto threshold=.95;//this is rather arbitrary
   
        auto filteredStrikes=filter(paddedStrikes, [&](const auto& x, const auto& i){
            return x<threshold;
        }, [&](const auto& x, const auto& i){
            return x;
        }, [&](const auto& x, const auto& i){
            return x;
        });

        auto filteredPrices=filter(optionPrices, [&](const auto& x, const auto& i){
            return paddedStrikes[i]<threshold;
        }, [&](const auto& x, const auto& i){
            return x-valOrZero(1-paddedStrikes[i]*discount);
        }, [&](const auto& x, const auto& i){
            return log(x);
        });
        auto sLow=spline::generateSpline(std::move(std::get<0>(filteredStrikes)), std::move(std::get<0>(filteredPrices)));
        auto sHigh=spline::generateSpline(std::move(std::get<1>(filteredStrikes)), std::move(std::get<1>(filteredPrices)));

        return [sLow=std::move(sLow), sHigh=std::move(sHigh), valOrZero=std::move(valOrZero), discount=std::move(discount), threshold=std::move(threshold)](const auto& k){
            if(k<threshold){
                return sLow(k);
            }
            else {
                return exp(sHigh(k))-valOrZero(1-k*discount);
            }
        };
        
    }


    template<typename Strike, typename MarketPrice, typename AssetPrice, typename R, typename T>
    auto generateFOEstimate(const std::vector<Strike>& strikes, const std::vector<MarketPrice>& options,  const AssetPrice& stock, const R& r, const T& t, const Strike& minStrike, const Strike& maxStrike){
        auto discount=exp(-r*t);
        auto s=getOptionSpline(strikes, options, stock, discount, minStrike, maxStrike);
        return [
            spline=std::move(s), 
            minStrike=transformPrice(minStrike, stock), 
            maxStrike=transformPrice(maxStrike, stock),
            discount=std::move(discount),
            t=std::move(t)
        ](int N, const auto& uArray){
            const auto xMin=log(discount*minStrike);
            const auto xMax=log(discount*maxStrike);
            auto valOrZero=[](const auto& v){
                return v>0.0?v:0.0;
            };
            return DFT(uArray, [&](const auto& u, const auto& x, const auto& index){
                const auto expX=exp(x);
                const auto strike=expX/discount;
                const auto offset=1.0-exp(x);
                const auto optionPrice=spline(strike);
                return valOrZero(optionPrice);
            }, [&](const auto& u, const auto& value){
                const auto front=u*cmp*(1.0+u*cmp);
                return log(1.0+front*value);
            }, xMin, xMax, N);
        };   

    }
    
    template<typename PhiHat, typename LogCfFN, typename DiscreteU>
    auto getObjFn_arr(PhiHat&& phiHattmp, LogCfFN&& cfFntmp, DiscreteU&& uArraytmp){
        return [
            phiHat=std::move(phiHattmp), 
            cfFn=std::move(cfFntmp), 
            uArray=std::move(uArraytmp)
        ](const auto& arrayParam){
            return futilities::sum(uArray, [&](const auto& u, const auto& index){
                return std::norm(
                    phiHat[index]-cfFn(std::complex<double>(1.0, u), arrayParam)
                )/(double)uArray.size();         
            });
        };
    }

}

#endif