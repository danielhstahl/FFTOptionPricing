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
    template<typename Index, typename Number>
    auto getUDomain(const Number& uMax, const Number& du, const Index& index){
        auto myU=du*index;
        return myU>uMax?myU-2.0*uMax:myU;
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
    auto FSTS(
        const Index& numSteps, 
        const Number& xmin, 
        const Number& xmax, 
        Discount&& discount, 
        Payoff&& payoff, 
        CF&& cf
    ){
        auto dx=computeDX(numSteps, xmin, xmax);
        //auto vmax=.5*(numSteps-1.0)/(xmax-xmin);
        auto vmax=M_PI/dx;
        auto du=2.0*vmax/numSteps;
        //auto vmin=0.0;
        auto vmin=du-vmax;
        //auto du=2.0*vmax/(numSteps-1.0);
        auto incrementedPayoff=ifft(
            futilities::for_each_parallel(
                fft(
                    futilities::for_each_parallel(0, numSteps, [&](const auto& index){
                        return payoff(getDomain(-xmax, dx, index))*exp(
                        std::complex<double>(0, -getDomain(0.0, dx*vmin, index)));
                    })
                ), 
                [&](const auto& val, const auto& index){
                    return val*cf(std::complex<double>(0, getDomain(vmin, du, index)));
                }
            )
        );
        return futilities::for_each_parallel(0, numSteps, [&](const auto& index){
            return discount(getDomain(-xmax, dx, index))*
            (incrementedPayoff[index]*exp(
                std::complex<double>(0, getDomain(0.0, vmin*dx, index)))
                ).real()/numSteps;
        });
    }


    template<typename Index, typename Number, typename Discount, typename Payoff, typename CF>
    auto FSTS(
        const Index& numSteps, 
        const Number& xmax, 
        Discount&& discount, 
        Payoff&& payoff, 
        CF&& cf
    ){
        return FSTS(numSteps, -xmax, xmax, discount, payoff, cf);
    }

}

#endif