| [Linux][lin-link] | [Windows][win-link] | [Codecov][cov-link] |
| :---------------: | :-----------------: | :-------------------: |
| ![lin-badge]      | ![win-badge]        | ![cov-badge]          |

[lin-badge]: https://travis-ci.org/phillyfan1138/FFTOptionPricing.svg?branch=master "Travis build status"
[lin-link]:  https://travis-ci.org/phillyfan1138/FFTOptionPricing "Travis build status"
[win-badge]: https://ci.appveyor.com/api/projects/status/i7agjioyxflo0xgq?svg=true "AppVeyor build status"
[win-link]:  https://ci.appveyor.com/project/phillyfan1138/fftoptionpricing "AppVeyor build status"
[cov-badge]: https://codecov.io/gh/phillyfan1138/FFTOptionPricing/branch/master/graph/badge.svg
[cov-link]:  https://codecov.io/gh/phillyfan1138/FFTOptionPricing



## FFT Option Pricing

This repository does Carr Madan, Fourier Space Time Stepping Algorithms, and Fang Oosterlee pricing from <a href="http://engineering.nyu.edu/files/jcfpub.pdf">Carr and Madan</a>, <a href="https://tspace.library.utoronto.ca/bitstream/1807/19300/1/Surkov_Vladimir_200911_PhD_Thesis.pdf">Sukov Vladimir</a>, and <a href="http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf">Fang Oosterlee</a>. 
Requires my <a href="https://github.com/phillyfan1138/FunctionalUtilities">Functional Utilities</a>, my <a href="https://github.com/phillyfan1138/CharacteristicFunctions">Characteristic Functions</a>, my <a href="https://github.com/phillyfan1138/FangOost">Fang Oosterlee </a>, and my <a href="https://github.com/phillyfan1138/RungeKutta">Runge Kutta</a> library.  

Benefits to FSTS: It only requires a CF and a payoff.  No further manipulations necessary.  Can also price American options.  Disadvantages: slower.  Roughly twice as slow as Carr and Madan for Europeans.  Prices in log asset instead of log strike; so not as easy to price all strikes for a given stock price. 

Benefits to Carr Madan: Ubiquitous. Prices in log strike.  Disadvantages: requires algebraic manipulations for each payoff.  Does not converge quickly.

Benefits to Fang Oosterlee:  Fast convergence in u.  U and X are separate and so can price only the options available in the market without interpolation.  Prices in either asset or strike.  Disadvantages: requires algebraic manipulations for each payoff.   For a given size of "X", tends to be slower than Carr Madan.  Benefit is that "X" can be whatever the market provides. 

## Examples

See the [tests](./test.cpp) for examples on how to use each of these algorithms.

## Generate charts and calibration documentation

For the calibration aspect, run `make generateCharts` and then run `./generateCharts`.  Then open `OptionCalibration.Rnw` in RStudio and compile as PDF.  