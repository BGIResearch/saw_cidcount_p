cc := g++
src := ./src
all=a.out
DFLAG=-lhdf5 -ldeflate
$(all):
	$(cc) -fopenmp -std=c++11 main.cpp -o $@ $(DFLAG)