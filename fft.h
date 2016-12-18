#ifndef __FFT_H_INCLUDED__
#define __FFT_H_INCLUDED__
#include <complex>
#include <iostream>
#include <vector>
#include <cmath>
#define _USE_MATH_DEFINES
typedef std::complex<double> Complex;
typedef std::vector<Complex> CArray;
 
//CArray fft(CArray&&);
//CArray ifft(CArray&&);
//void fft(int, std::vector<Complex>&);
#include "fft.hpp"
#endif
