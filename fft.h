#ifndef __FFT_H_INCLUDED__
#define __FFT_H_INCLUDED__
#include <complex>
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

typedef std::complex<double> Complex;
typedef std::vector<Complex> CArray;

#include "fft.hpp"
#endif
