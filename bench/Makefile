BENCH_PRGS=bench_elimination 

CPUCYCLES_DIR=./cpucycles-20060326

all: $(BENCH_PRGS)


bench_elimination: cpucycles.o	bench_elimination.cc
	g++ -O2 -DNDEBUG -I$(SAGE_ROOT)/local/include  -L/$(SAGE_ROOT)/local/lib -I../ -lm4rie -lm4ri -lgivaro -lntl -lgmpxx -lgmp -lm -lstdc++ cpucycles.o bench_elimination.cc -o $@

clean:
	rm -f $(BENCH_PRGS)
	rm -f *.o cpucycles.h
	rm -f $(CPUCYCLES_DIR)/cpucycles.o
	rm -f $(CPUCYCLES_DIR)/cpucycles.h



cpucycles.h: cpucycles.o

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; cp cpucycles.o ..; cp cpucycles.h ..; cd ..)
