#ifndef OPTIONPRICING_H
#define OPTOINPRICING_H
#include <complex>
#include <vector>
#include "FunctionalUtilities.h"
#include "CharacteristicFunctions.h"
#include "RungeKutta.h"
#include "fft.h"
namespace optionprice{
    template<typename Index, typename Number>
    auto computeDX(const Index& xDiscrete, const Number& xMin,const Number& xMax){
        return (xMax-xMin)/(xDiscrete-1.0);
    }
    template<typename Index, typename Number>
    auto getDomain(const Number& xMin, const Number& dx, const Index& index){
        return xMin+index*dx;
    }
    template<typename Number1, typename Number2>
    auto getRealDomain(const Number1& startPoint, const Number2& logDomain){
        return startPoint*exp(logDomain);
    }

    template<typename Number>
    std::vector<Number> computeXRange(int xDiscrete, const Number& xMin, const Number& xMax){
        return futilities::for_emplace_back(xMin, xMax, xDiscrete, [](const auto& val){
            return val;
        });
    }



    /**
        Fourier Space-Time Stepping algorithm
        returns in log domain
    */

    template<typename Index, typename Number, typename Discount, typename Payoff, typename CF>
    auto FSTS(const Index& numSteps, const Number& xmin, const Number& xmax, /*const Number& Maturity,*/ Discount&& discount, Payoff&& payoff, CF&& cf){
        auto dx=computeDX(numSteps, xmin, xmax);
        auto vmax=M_PI/dx;
        auto du=2.0*vmax/numSteps;
        auto vmin=du-vmax;
        auto incrementedPayoff=ifft(
            futilities::for_each_parallel(
                fft(
                    futilities::for_each_parallel(0, numSteps, [&](const auto& index){
                        return payoff(getDomain(xmin, dx, index))*exp(std::complex<double>(0, -getDomain(0.0, vmin*dx, index)));
                    })
                ), 
                [&](const auto& val, const auto& index){
                    auto u=std::complex<double>(getDomain(vmin, du, index), 0);
                    return val*cf(u);
                }
            )
        );
        //std::vector<priceAndUnderlying> result(numSteps);
        return futilities::for_each_parallel(0, numSteps, [&](const auto& index){
            return discount(getDomain(xmin, dx, index))*(incrementedPayoff[index]*exp(std::complex<double>(0, getDomain(0.0, vmin*dx, index)))).real()/numSteps;
        });
    }

}

#endif