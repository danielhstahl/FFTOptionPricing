INCLUDES=-I ../FunctionalUtilities -I ../RungeKutta -I ../CharacteristicFunctions
test:test.o
	g++ -std=c++14 -O3  test.o $(INCLUDES) -o test -fopenmp
test.o: main.cpp fft.h payoffs.h fft.hpp FSTS.hpp payoffs.hpp
	g++ -std=c++14 -O3  -c test.cpp  $(INCLUDES) -fopenmp
clean:
	-rm *.o test
