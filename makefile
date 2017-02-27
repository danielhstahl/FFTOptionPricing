INCLUDES=-I ../FunctionalUtilities -I ../RungeKutta -I ../CharacteristicFunctions -I ../FangOost
test:test.o
	g++ -std=c++14 -O3  -pthread --coverage test.o $(INCLUDES) -o test -fopenmp
test.o: test.cpp fft.h fft.hpp OptionPricing.h
	g++ -std=c++14 -O3  -pthread --coverage -c test.cpp  $(INCLUDES) -fopenmp 
clean:
	-rm *.o test
