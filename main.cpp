#include "Complex.h"
#include <cmath>
#include <iostream>
#include "FSTS.h"
#include "payoffs.h"
#include "CF.h"
#include "fft.h"
#include <vector>
#include "CarrMadan.h"
int main(){
  double r=.03;
  double sigma=.3;
  double t=1;
  double S0=40;
  double K=40;
  double discount=exp(-r*t);
  double drift=r-sigma*sigma*.5;
  auto cf=[&](Complex& u, double t){
    return BMWithDrift(u, drift, sigma, t);
  };
  int numSteps=pow(2, 8);
  auto payoff=[&](double y){
    return Call(S0, y, K);
  };
  auto payoffCF=[&](double v, double t, double alpha, auto& cf){
    return Call(v, t, alpha, cf);
  };
  double xmin=-5;
  double xmax=5;
  FSTS fsts(numSteps, xmin, xmax);
  std::vector<priceAndUnderlying> price=fsts.OptionPrice(t, discount, payoff, cf);
  int cutoff=.6*numSteps;
  for(int i=0; i<cutoff; ++i){
    //std::cout<<S0*exp(price[i].underlying)<<", "<<price[i].price<<std::endl;
  }

  CarrMadan carrm(numSteps);
  std::vector<priceAndStrike> priceC=carrm.OptionPrice(t, discount, payoffCF, cf);

  for(int i=0; i<cutoff; ++i){
    std::cout<<S0*exp(priceC[i].strike)<<", "<<priceC[i].price*S0<<std::endl;
  }


}
