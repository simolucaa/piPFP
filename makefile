# compilation flags
CXX_FLAGS=-std=c++20 -O3 -Wall -Wextra -g -march=native
CFLAGS=-O3 -Wall -std=c11 -g
CC=gcc

# main executables 
EXECS=newscan.x pscan.x
# executables not using threads (and therefore not needing the thread library)
EXECS_NT=newscanNT.x 

# targets not producing a file declared phony
.PHONY: all clean lcp

all: $(EXECS) $(EXECS_NT) newscanNT.x newscan.x newscan_faster.x newscan_fasterNT.x pscan.x newscan_faster_growthNT.x newscan_faster_growth.x newscan_faster_growth2.x newscan_faster_growth2NT.x


newscanNT.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS

newscan.x: newscan.cpp newscan.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan.cpp malloc_count.o utils.o xerrors.o -ldl -pthread

newscan_faster.x: newscan_faster.cpp newscan_faster.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan_faster.cpp malloc_count.o utils.o xerrors.o -ldl -pthread -ltbb

newscan_fasterNT.x: newscan_faster.cpp newscan_faster.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan_faster.cpp malloc_count.o utils.o xerrors.o -ldl -pthread -ltbb -DNOTHREADS

newscan_faster_growth.x: newscan_faster_growth.cpp newscan_faster_growth.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan_faster_growth.cpp malloc_count.o utils.o xerrors.o -ldl -pthread -ltbb

newscan_faster_growthNT.x: newscan_faster_growth.cpp newscan_faster_growth.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan_faster_growth.cpp malloc_count.o utils.o xerrors.o -ldl -pthread -ltbb -DNOTHREADS

newscan_faster_growth2.x: newscan_faster_growth2.cpp newscan_faster_growth2.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan_faster_growth2.cpp malloc_count.o utils.o xerrors.o -ldl -pthread -ltbb

newscan_faster_growth2NT.x: newscan_faster_growth2.cpp newscan_faster_growth2.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan_faster_growth2.cpp malloc_count.o utils.o xerrors.o -ldl -pthread -ltbb -DNOTHREADS

pscan.x: pscan.cpp pscan.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ pscan.cpp malloc_count.o utils.o xerrors.o -ldl -pthread

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.x
