#ifndef OPTIONPRICING_H
#define OPTIONPRICING_H
#include <complex>
#include <vector>
#include <algorithm>
#include "FunctionalUtilities.h"
#include "CharacteristicFunctions.h"
#include "FangOost.h"
#include "RungeKutta.h"
#include "fft.h"
namespace optionprice{
    /**
    Used by FSTS method
    @xDiscrete Number of nodes used to discretize X
    @xMin minimum x
    @xMax maximum x
    @returns the distance between nodes
    */
    /*template<typename Index, typename Number>
    auto computeDX(const Index& xDiscrete, const Number& xMin,const Number& xMax){
        return (xMax-xMin)/(xDiscrete-1.0);
    }*/
    /**
    Used by FSTS method
    @xMin minimum x
    @xMax distance between nodes, eg from computeDX
    @index index of the node
    @returns value at the node index
    */
    /*template<typename Index, typename Number>
    auto getDomain(const Number& xMin, const Number& dx, const Index& index){
        return xMin+index*dx;
    }*/

    /**
    Used by FSTS method.  This wraps back around when reaches a certain level.  
    @uMax maximum u
    @du distance between nodes
    @index index of the node
    @returns value at the node index
    */
    template<typename Index, typename Number>
    auto getUDomain(const Number& uMax, const Number& du, const Index& index){
        auto myU=du*index;
        return myU>uMax?myU-2.0*uMax:myU;
    }

    /**
    Used by Carr Madan
    @ada the distance between nodes in the U domain
    @returns the maximum value of the U domain
    */
    template<typename Number>
    auto getMaxK(const Number& ada){
        return M_PI/ada;
    }
    /**
    Used by Carr Madan
    @numSteps the number of steps used to discretize the x and u domains
    @returns the distance between nodes in the X domain
    */
    template<typename Index, typename Number>
    auto getLambda(const Index& numSteps, const Number& b){
        return 2.0*b/numSteps;
    }


    /**
    Used by FSTS
    @xMin minimum X
    @xMax maximum X
    @K Strike price
    @numX number of nodes
    @returns vector of stock or asset prices.  Note that this is different from Carr Madan which prices in terms of a vector of strike prices
    */
    /*template<typename Index, typename Number>
    auto getFSTSUnderlying(const Number& xMin, const Number& xMax, const Number& K, const Index& numX){
        auto dx=fangoost::computeDX(numX, xMin, xMax);
        return futilities::for_each_parallel(0, numX, [&](const auto& index){
            return K*exp(getDomain(xMin, dx, index));
        });
    }*/
    /**
    Helper function for Fang Oosterlee
    @xMin minimum X
    @xMax maximum X
    @K Strike price
    @numX number of nodes
    @returns vector of stock or asset prices.  Note that this is different from Carr Madan which prices in terms of a vector of strike prices.  However, this can also price per strike; see getFangOostStrike
    */
    template<typename Index, typename Number>
    auto getStrikeUnderlying(const Number& xMin, const Number& xMax, const Number& K, const Index& numX){
        auto dx=fangoost::computeDX(numX, xMin, xMax);
        return futilities::for_each_parallel(0, numX, [&](const auto& index){
            return K*exp(fangoost::getX(xMin, dx, index));
        });
    }

    template<typename Index, typename Number>
    auto getFangOostKAtIndex(const Number& xMin, const Number& dK, const Number& S0, const Index& index){
        return S0*exp(-fangoost::getX(xMin, dK, index));
    }
    /**helper function for computing strike vector for 
    testing purposes.  Converts x back to strike*/
    template<typename Index, typename Number>
    auto getFangOostStrike(const Number& xMin, const Number& xMax, const Number& S0, const Index& numX){
        auto dx=fangoost::computeDX(numX, xMin, xMax);
        return futilities::for_each_parallel(0, numX, [&](const auto& index){
            //the exponent is negative because x=log(S0/k)->k=s0*exp(-x)
            return getFangOostKAtIndex(xMin, dx, S0, index);
        });
    }
    /**This function takes strikes and converts them
    into a vector in the x domain.  Intriguinely, I 
    don't have to sort the result...*/
    template< typename Array, typename Number>
    auto getXFromK(const Number& S0, const Array& K){
        return futilities::for_each_parallel_copy(
            K, 
            [&](const auto& val, const auto& index){
                return log(S0/val);
            }
        );
    }


    template<typename Index, typename Number>
    auto getCarrMadanKAtIndex(const Number& b, const Number& lambda, const Number& S0, const Index& index){
        return S0*exp(fangoost::getX(-b, lambda, index));
    }
    /**
    Used by CarrMadan
    @ada distance between nodes in the U domain
    @S0 stock or asset price
    @numX number of nodes
    @returns vector of strikes  Note that this is different from FSTS and FangOosterlee which prices in terms of a vector of stock or asset prices
    */
    template<typename Index, typename Number>
    auto getCarrMadanStrikes(const Number& ada, const Number& S0, const Index& numX){
        auto b=getMaxK(ada);
        auto lambda=getLambda(numX, b);
        return futilities::for_each_parallel(0, numX, [&](const auto& index){
            return getCarrMadanKAtIndex(b, lambda, S0, index);
        });
    }


    /**
        Fourier Space-Time Stepping algorithm 
        returns in log domain
        respect to S, not K
        https://tspace.library.utoronto.ca/bitstream/1807/19300/1/Surkov_Vladimir_200911_PhD_Thesis.pdf
        @numSteps Discrete steps in X domain
        @xmax maximum X (in log asset space around the strike)
        @discount function which takes the X values as input.  Can be as simple as e^-rt but can be useful for interest rate options
        @payoff payoff function which takes the log result.  
        @CF characteristic function of log x around the strike
        @returns vector of prices corresponding with the assets given in getFSTSUnderlying
    */
    template<typename Index, typename Number, typename Discount, typename Payoff, typename CF>
    auto FSTS(
        const Index& numSteps, 
        const Number& xmax, 
        Discount&& discount, 
        Payoff&& payoff, 
        CF&& cf
    ){
        auto dx=fangoost::computeDX(numSteps, -xmax, xmax);
        auto vmax=M_PI/dx;
        auto du=2.0*vmax/numSteps;
        auto incrementedPayoff=ifft(
            futilities::for_each_parallel(
                fft(
                    futilities::for_each_parallel(0, numSteps, [&](const auto& index){
                        return std::complex<double>(payoff(fangoost::getX(-xmax, dx, index)), 0);
                    })
                ), 
                [&](const auto& val, const auto& index){
                    return val*cf(std::complex<double>(0, getUDomain(vmax, du, index)));
                }
            )
        );
        return futilities::for_each_parallel(0, numSteps, [&](const auto& index){
            return discount(fangoost::getX(-xmax, dx, index))*incrementedPayoff[index].real()/numSteps;
        });
    }
    

    /**For Fang Oost (defined in the paper)*/
    template<typename A, typename B, typename C, typename D, typename U>
    auto chiK(const A& a, const B& b, const C& c, const D& d, const U& u){
        auto iterS=[&](const auto& x){
            return u*(x-a);
        };
        auto expD=exp(d);
        auto expC=exp(c);
        return (cos(iterS(d))*expD-cos(iterS(c))*expC+u*sin(iterS(d))*expD-u*sin(iterS(c))*expC)/(1.0+u*u);
    }
    /**For Fang Oost (defined in the paper)*/
    template<typename A, typename B, typename C, typename D, typename U, typename Index>
    auto phiK(const A& a, const B& b, const C& c, const D& d, const U& u, const Index& k){
        auto iterS=[&](const auto& x){
            return u*(x-a);
        };
        return k==0?d-c:(sin(iterS(d))-sin(iterS(c)))/u;
    }
        
    /**
        Fang Oosterlee Approach for an option using Put as the main payoff (better accuracy than a call...use put call parity to get back put).
        Note that Fang Oosterlee's approach works well for a smaller 
        of discrete strike prices such as those in the market.  The 
        constraint is that the smallest and largest values in the x domain
        must be relatively far from the middle values.  This can be 
        "simulated" by adding small and large "K" synthetically.  Due to
        the fact that Fang Oosterlee is able to handle this well, the 
        function takes a vector of strike prices with no requirement that
        the strike prices be equidistant.  All that is required is that
        they are sorted largest to smallest.

        returns in log domain
        http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf
        @xValues x values derived from strikes
        @numUSteps number of steps in the complex domain (independent 
        of number of x steps)
        @mOutput a function which determines whether the output is a 
        call or a put.  
        @CF characteristic function of log x around the strike
        @returns vector of prices corresponding with the strikes 
        provided by FangOostCall or FangOostPut
    */
    template<typename Index, typename Array,  typename CF, typename OutputManipulator>
    auto FangOostGeneric(
        Array&& xValues,
        const Index& numUSteps, 
        OutputManipulator&& mOutput,
        CF&& cf
     ){
         //x goes from log(S0/kmin) to log(S0/kmax)
        auto xMin=xValues.front();
        auto xMax=xValues.back();
        return futilities::for_each_parallel(
            fangoost::computeExpectationVectorLevy(
                xValues, numUSteps, cf, 
                [&](const auto& u, const auto& x, const auto& k){
                    return phiK(xMin, xMax, xMin, 0.0, u, k)-chiK(xMin, xMax, xMin, 0.0, u);//used for put
                }
            ), 
            mOutput
        );
    }

    

    /**
        Fang Oosterlee Approach for a PUT (better accuracy than a call...use put call parity to get back put)
        http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf
        @S0 value of stock
        @K stl container containing strikes to be priced..must be in DESCENDING order
        @numUSteps number of steps in complex domain
        @discount constant which is used for discounting
        @CF characteristic function of log x around the strike
        @returns vector of put prices corresponding with the K
    */
    template<typename Index, typename Array, typename Number,  typename CF>
    auto FangOostPut(
        const Number& S0,
        const Array& K,
        const Index& numUSteps,     
        const Number& discount,
        CF&& cf
     ){
        return FangOostGeneric(
            getXFromK(S0, K),
            numUSteps,  
            [&](const auto& val, const auto& index){
                return val*discount*K[index];
            }, 
            cf
        );
    }
    /**
        Fang Oosterlee Approach for a call 
        http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf
        @S0 value of stock
        @K stl container containing strikes to be priced..must be in DESCENDING order
        @numUSteps number of steps in complex domain
        @discount constant which is used for discounting
        @CF characteristic function of log x around the strike
        @returns vector of call prices corresponding with the K
    */
    template<typename Index, typename Array, typename Number, typename CF>
    auto FangOostCall(
        const Number& S0,
        const Array& K,
        const Index& numUSteps,         
        const Number& discount,
        CF&& cf
     ){
        return FangOostGeneric(
            getXFromK(S0, K),
            numUSteps, 
            [&](const auto& val, const auto& index){
                return (val-1.0)*discount*K[index]+S0;
            }, 
            cf
        );
    }
    
    /**
    Used for Carr-Madan call option http://engineering.nyu.edu/files/jcfpub.pdf
    */
    template<typename Number1, typename Number2,  typename CF>
    auto CallAug(const Number1& v, const Number2& alpha, CF& cf){ //used for Carr-Madan approach...v is typically complex
        return cf(v+(alpha+1.0))/(alpha*alpha+alpha+v*v+(2*alpha+1)*v);
    }
    
    /**
        Carr Madan Approach for using a CALL as the fundamental asset
        respect to K not S
        http://engineering.nyu.edu/files/jcfpub.pdf
        @numSteps Discrete steps in X domain
        @ada distance between nodes in U domain
        @alpha a variable used to keep integrand from 0
        @S0 asset value
        @discount constant which is used for discounting
        @mOutput function of fundamental asset.  For a call, this just returns the value
        @CF characteristic function of log x around the strike
        @returns vector of prices corresponding with the strikes given in getCarrMadanStrikes
    */
    template<typename Index, typename Number,  typename CF, typename OutputManipulator>
    auto CarrMadanGeneric(
        const Index& numSteps, 
        const Number& ada, 
        const Number& alpha, 
        const Number& S0, 
        const Number& discount, 
        OutputManipulator&& mOutput,
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
            return mOutput(S0*cmplVector[index].real()*exp(-alpha*fangoost::getX(-b, lambda, index))*ada/(M_PI*3.0), index);
        });
    }

    template<typename Index, typename Number,  typename CF>
    auto CarrMadanCall(
        const Index& numSteps, 
        const Number& ada, 
        const Number& alpha, 
        const Number& S0, 
        const Number& discount, 
        CF&& cf
    ){
        return CarrMadanGeneric(numSteps, ada, alpha, S0, discount, 
        [&](const auto& val, const auto& index){
            return val;
        }, [&](const auto& v, const auto& alpha){
            return optionprice::CallAug(v, alpha, cf);
        });
        
    }
    template<typename Index, typename Number,  typename CF>
    auto CarrMadanPut(
        const Index& numSteps, 
        const Number& ada, 
        const Number& alpha, 
        const Number& S0, 
        const Number& discount, 
        CF&& cf
    ){
        auto b=getMaxK(ada);
        auto lambda=getLambda(numSteps, b);
        return CarrMadanGeneric(numSteps, ada, alpha, S0, discount, 
        [&](const auto& val, const auto& index){
            return val-S0+getCarrMadanKAtIndex(b, lambda, S0, index)*discount;
        }, [&](const auto& v, const auto& alpha){
            return optionprice::CallAug(v, alpha, cf);
        });
    }
    

    /**
        Carr Madan Approach for a CALL 
        respect to K not S
        http://engineering.nyu.edu/files/jcfpub.pdf
        @numSteps Discrete steps in X domain
        @ada distance between nodes in U domain
        @S0 asset value
        @discount constant which is used for discounting
        @CF characteristic function of log x around the strike
        @returns vector of prices corresponding with the strikes given in getCarrMadanStrikes
    */
    template<typename Index, typename Number,  typename CF>
    auto CarrMadanCall(
        const Index& numSteps, 
        const Number& ada, 
        const Number& S0,
        const Number& discount, 
        CF&& augCF
    ){
        auto alpha=1.5;
        return CarrMadanCall(numSteps, ada, alpha, S0, discount, augCF);
    }
    /**
        Carr Madan Approach for a Put 
        respect to K not S
        http://engineering.nyu.edu/files/jcfpub.pdf
        @numSteps Discrete steps in X domain
        @ada distance between nodes in U domain
        @S0 asset value
        @discount constant which is used for discounting
        @CF characteristic function of log x around the strike
        @returns vector of prices corresponding with the strikes given in getCarrMadanStrikes
    */
    template<typename Index, typename Number,  typename CF>
    auto CarrMadanPut(
        const Index& numSteps, 
        const Number& ada, 
        const Number& S0,
        const Number& discount, 
        CF&& augCF
    ){
        auto alpha=1.5;
        return CarrMadanPut(numSteps, ada, alpha, S0, discount, augCF);
    }
}

#endif
