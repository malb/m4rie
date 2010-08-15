BENCH_PRGS=bench_elimination bench_multiplication

CPUCYCLES_DIR=./cpucycles-20060326

all: $(BENCH_PRGS)

CXXFLAGS=-I$(SAGE_ROOT)/local/include  -I../ -O2 -DNDEBUG 
LDFLAGS=-L/$(SAGE_ROOT)/local/lib -I../ -lm4rie -lm4ri -lgivaro -lntl -lgmpxx -lgmp -lm -lstdc++

bench_elimination: cpucycles.o	bench_elimination.cc
	g++ $(CXXFLAGS) $(LDFLAGS) cpucycles.o bench_elimination.cc -o $@

bench_multiplication: cpucycles.o bench_multiplication.cc
	g++ $(CXXFLAGS) $(LDFLAGS) cpucycles.o bench_multiplication.cc -o $@

clean:
	rm -f $(BENCH_PRGS)
	rm -f *.o cpucycles.h
	rm -f $(CPUCYCLES_DIR)/cpucycles.o
	rm -f $(CPUCYCLES_DIR)/cpucycles.h
	rm -rf callgrind.*



cpucycles.h: cpucycles.o

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; cp cpucycles.o ..; cp cpucycles.h ..; cd ..)
