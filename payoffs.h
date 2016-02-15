#ifndef __PAYOFFS_H_INCLUDED__
#define __PAYOFFS_H_INCLUDED__
#include <cmath>
#include "Complex.h"
double Call(
  double,//initial price
  double,//result of log
  double//strike
);
double Put(
  double,//initial price
  double,//result of log
  double//strike
);
Complex Call(
  double, //v
  double, //t
  double, //alpha
  auto& //cf
);
#include "payoffs.hpp"
#endif
