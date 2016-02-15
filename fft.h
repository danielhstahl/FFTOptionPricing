#ifndef __FFT_H_INCLUDED__
#define __FFT_H_INCLUDED__
#include "Complex.h"
#include <iostream>
#include <vector>
#include <cmath>
#define _USE_MATH_DEFINES
void fft(std::vector<Complex>&);
void ifft(std::vector<Complex>&);
//void fft(int, std::vector<Complex>&);
#include "fft.hpp"
#endif
