# compilation flags
CXX_FLAGS=-std=c++20 -O3 -Wall -Wextra -g -march=native
CXX_LIBS=-ldl -pthread -I/opt/homebrew/Cellar/tbb/2022.0.0/include -L/opt/homebrew/Cellar/tbb/2022.0.0/lib -ltbb
CFLAGS=-O3 -Wall -std=c11 -g
CC=gcc
CCX=g++

# targets not producing a file declared phony
.PHONY: all clean 

all: piPFP.x piPFP_NT.x piPFP_growth.x piPFP_growth_NT.x


piPFP.x: piPFP.cpp piPFP.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ piPFP.cpp malloc_count.o utils.o xerrors.o $(CXX_LIBS)

piPFP_NT.x: piPFP.cpp piPFP.hpp malloc_count.o utils.o xerrors.o
	$(CXX) $(CXX_FLAGS) -o $@ piPFP.cpp malloc_count.o utils.o xerrors.o $(CXX_LIBS) -DNOTHREADS

piPFP_growth.x: piPFP_growth.cpp piPFP_growth.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ piPFP_growth.cpp malloc_count.o utils.o xerrors.o $(CXX_LIBS)

piPFP_growth_NT.x: piPFP_growth.cpp piPFP_growth.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ piPFP_growth.cpp malloc_count.o utils.o xerrors.o $(CXX_LIBS) -DNOTHREADS

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.x
