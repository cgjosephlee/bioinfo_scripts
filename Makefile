CXX=g++-7
CXXFLAGS=-g -Wall -O2 -std=c++17
LIBS=-lz
PROG=fastq_stats

all:$(PROG)

fastq_stats: fastq_stats.cpp
	$(CXX) $(CXXFLAGS) -I./lib -o $@ $< $(LIBS)
