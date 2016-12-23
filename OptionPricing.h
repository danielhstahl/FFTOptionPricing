#ifndef OPTIONPRICING_H
#define OPTOINPRICING_H
#include <complex>
#include <vector>
#include "FunctionalUtilities.h"
#include "CharacteristicFunctions.h"
#include "FangOost.h"
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

    template<typename Index, typename Number>
    auto computeXRange(const Index& xDiscrete, const Number& xMin, const Number& xMax){
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
    Fang Oosterlee approach
    
    */
    template<typename Number, typename Index>
    auto chiK(const Number& a, const Number& b, const Number& c, const Number& d, const Index& k){
        auto efficientSquare=[](const auto& num){
            return num*num;
        };
        auto kPiBA=k*M_PI/(b-a);
        auto iterS=[&](const auto& x){
            return kPiBA*(x-a);
        }
        auto coef=1.0/(1+efficientSquare(kPiBA));
        auto expD=exp(d);
        auto expC=exp(c);
        return coef*(cos(iterS(d))*expD-cos(iterS(c))*expC+kPiBA*sin(iterS(d)*expD-kPiBA*sin(iterS(c))*expC);
    }
    template<typename Number, typename Index>
    auto phiK(const Number& a, const Number& b, const Number& c, const Number& d, const Index& k){
   
        auto kPiBA=k*M_PI/(b-a);
        auto iterS=[&](const auto& x){
            return kPiBA*(x-a);
        }
        return k==0?d-c:(sin(iterS(d))-sin(iterS(c)))/kPiBA);
    }
    
    /**
    For a call, the "VK" is (2/(b-a))*K(chiK(0, b)-phiK(0, b))
    */
    template<typename Index, typename Number,  typename CF>
    auto FangOost(
        const Index& numXSteps, 
        const Index& numUSteps, 
        const Number& xMin, 
        const Number& xMax, 
        const Number& K, 
        const Number& discount
        CF&& cf
     ){
    
        //int xDiscrete, int uDiscrete,  const Number& xMin, const Number& xMax, auto&& fnInv, auto&& vK
        return fangoost::computeInv(numXSteps, numUSteps, xMin, xMax, cf, 
           /* [&](const auto& u, const auto& x){
                return K*(chiK(xMin, xMax, 0.0, xMax, index)-phiK(xMin, xMax, 0.0, xMax, index))*2.0/(xMax-xMin);
            });*/
        
        
        futilities::for_each_parallel(0, numSteps, [&](const auto& index){
            return K*(chiK(xMin, xMax, 0.0, xMax, index)-phiK(xMin, xMax, 0.0, xMax, index))*2.0/(xMax-xMin);
        }));
        
    }
    template<typename Index, typename Number,  typename CF>
    auto FangOost(
        const Index& numXSteps, 
        const Index& numUSteps, 
        const Number& xMax, 
        const Number& K, 
        const Number& discount
        CF&& cf
     ){
        return FangOost(numXSteps, numUSteps, -xMax, xMax, K, discount, cf);
     }
    
    /**
    Used for Carr-Madan call option http://engineering.nyu.edu/files/jcfpub.pdf
    */
    template<typename Number1, typename Number2,  typename CF>
    auto CallAug(const Number1& v, const Number2& alpha, CF& cf){ //used for Carr-Madan approach...v is typically complex
        //auto u=;
        return cf(v+(alpha+1.0))/(alpha*alpha+alpha+v*v+(2*alpha+1)*v);
    }
    template<typename Number>
    auto getMaxK(const Number& ada){
        return M_PI/ada;
    }
    template<typename Index, typename Number>
    auto getLambda(const Index& numSteps, const Number& b){
        return 2.0*b/numSteps;
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
        auto lambda=getLambda(numSteps, b);
        auto adjustFirst=[](auto&& array){array[0]=array[0]/2.0;return std::move(array);};
        auto cmplVector=fft(
            adjustFirst(
                futilities::for_each_parallel(0, numSteps, [&](const auto& index){
                    auto pm=index%2==0?-1.0:1.0;
                    return discount*augCF(std::complex<double>(0, index*ada), alpha)*exp(std::complex<double>(0, b*index*ada))*(3.0+pm);
                })
            )
        );
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
