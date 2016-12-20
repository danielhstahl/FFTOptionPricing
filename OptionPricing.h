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
        respect to S, not K
        https://tspace.library.utoronto.ca/bitstream/1807/19300/1/Surkov_Vladimir_200911_PhD_Thesis.pdf
    */

    template<typename Index, typename Number, typename Discount, typename Payoff, typename CF>
    auto FSTS(
        const Index& numSteps, 
        const Number& xmax, 
        Discount&& discount, 
        Payoff&& payoff, 
        CF&& cf
    ){
        auto dx=computeDX(numSteps, -xmax, xmax);
        auto vmax=M_PI/dx;
        auto du=2.0*vmax/numSteps;
        auto incrementedPayoff=ifft(
            futilities::for_each_parallel(
                fft(
                    futilities::for_each_parallel(0, numSteps, [&](const auto& index){
                        return std::complex<double>(payoff(getDomain(-xmax, dx, index)), 0);
                    })
                ), 
                [&](const auto& val, const auto& index){
                    return val*cf(std::complex<double>(0, getUDomain(vmax, du, index)));
                }
            )
        );
        return futilities::for_each_parallel(0, numSteps, [&](const auto& index){
            return discount(getDomain(-xmax, dx, index))*incrementedPayoff[index].real()/numSteps;
        });
    }


    /**
    Used for Carr-Madan call option http://engineering.nyu.edu/files/jcfpub.pdf
    */
    template<typename Number1, typename Number2,  typename CF>
    auto CallAug(const Number1& v, const Number2& alpha, CF& cf){ //used for Carr-Madan approach...v is typically complex
        auto u=v+(alpha+1.0);
        return cf(u)/(alpha*alpha+alpha+v*v+(2*alpha+1)*v);
    }
    template<typename Number>
    auto getMaxK(const Number& ada){
        return M_PI/ada;
    }
    /**
    Carr-Madan: with respsect to K, not S

    */
    template<typename Index, typename Number,  typename CF>
    auto CarrMadan(
        const Index& numSteps, 
        const Number& ada, 
        const Number& alpha, 
        const Number& discount, 
        CF&& augCF
    ){
        auto b=getMaxK(ada);
        auto lambda=2.0*b/numSteps;
        auto adjustFirst=[](auto&& array){
            array[0]=array[0]/2.0;
            return std::move(array);
        };
        auto cmplVector=fft(adjustFirst(futilities::for_each_parallel(0, numSteps, [&](const auto& index){
            auto pm=index%2==0?-1.0:1.0;
            return discount*augCF(std::complex<double>(0, index*ada), alpha)*exp(std::complex<double>(0, b*index*ada))*(3.0+pm);
        })));
        return futilities::for_each_parallel(0, numSteps, [&](const auto& index){
            return cmplVector[index].real()*exp(-alpha*getDomain(-b, lambda, index))*ada/(M_PI*3.0);
        });
    }
    template<typename Index, typename Number,  typename CF>
    auto CarrMadan(
        const Index& numSteps, 
        const Number& ada, 
        const Number& discount, 
        CF&& augCF
    ){
        auto alpha=1.5;
        return CarrMadan(numSteps, ada, alpha, discount, augCF);
    }
}

#endif