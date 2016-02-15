#ifndef __PAYOFFS_H_INCLUDED__
#define __PAYOFFS_H_INCLUDED__
#include <cmath>
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
#include "payoffs.hpp"
#endif
