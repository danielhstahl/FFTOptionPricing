#ifndef __CARRMADAN_H_INCLUDED__
#define __CARRMADAN_H_INCLUDED__
#include "fft.h"
#include "Complex.h"
#include <iostream>
struct priceAndStrike{
  double price;
  double strike;
};
class CarrMadan{ //carr-madan option pricing
private:
  int numSteps;
  double ada;
  double lambda;
  double b;
  double alpha;
  void setValues(int, double, double);
public:
  CarrMadan(int, double, double);//numSTeps, ada, alpha
  CarrMadan(int);//numSteps
  std::vector<priceAndStrike> OptionPrice(
    double,//t
    double,//discount
    auto&, //augCF
    auto& //characteristic function
  );
};

#include "CarrMadan.hpp"
#endif
