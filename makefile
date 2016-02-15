LDFLAGS=-L ../Complex -lComplex
INCLUDES=-I ../Complex
optionPricing:main.o
	g++ -std=c++14 -O3  main.o $(LDFLAGS) $(INCLUDES) -o optionPricing -fopenmp
main.o: main.cpp fft.h FSTS.h payoffs.h CF.h fft.hpp FSTS.hpp payoffs.hpp CF.hpp CarrMadan.h CarrMadan.hpp
	g++ -std=c++14 -O3  -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	-rm *.o optionPricing
