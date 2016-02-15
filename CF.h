#ifndef __CF_H_INCLUDED__
#define __CF_H_INCLUDED__
#include "Complex.h"
#include <cmath>

Complex BMWithDrift(
  Complex&, //u
  double, //drift (eg, r-sigma*sigma*.5)
  double, //volatility (eg, sigma)
  double //t
);
#include "CF.hpp"
#endif
