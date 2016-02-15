#ifndef __FSTS_H_INCLUDED__
#define __FSTS_H_INCLUDED__
#include "fft.h"
#include <iostream>
struct priceAndUnderlying{
  double price;
  double underlying;
};
class FSTS{ //fourier space time steps
private:
  int numSteps;
  double xmin;
  double xmax;
  double dx;
  double du;
  double vmax;
public:
  FSTS(int, double, double);//numSTeps, xmin, xmax
  std::vector<priceAndUnderlying> OptionPrice(
    double,//t
    double,//discount
    auto&, //payoff function
    auto& //characteristic function
  );
};

#include "FSTS.hpp"
#endif
