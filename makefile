INCLUDES=-I ../FunctionalUtilities -I ../RungeKutta -I ../CharacteristicFunctions -I ../FangOost -I ../MonteCarlo  -I../cuckoo_search 


GCCVAL=g++


UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	GCCVAL=g++-7
endif

test:test.o
	$(GCCVAL) -std=c++14 -O3  -pthread --coverage test.o $(INCLUDES) -o test -fopenmp
test.o: test.cpp fft.h fft.hpp OptionPricing.h OptionCalibration.h monotone_spline.h
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage -c test.cpp  $(INCLUDES) -fopenmp

generateCharts:generateCharts.o
	$(GCCVAL) -std=c++14 -O3  -pthread --coverage generateCharts.o $(INCLUDES) -o generateCharts -fopenmp
generateCharts.o: generateCharts.cpp fft.h fft.hpp OptionPricing.h OptionCalibration.h monotone_spline.h
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage -c generateCharts.cpp  $(INCLUDES) -fopenmp

clean:
	-rm *.o test *.csv *.pdf


