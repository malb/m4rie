CPUCYCLES_DIR=./cpucycles-20060326
MYLDFLAGS=-L../.libs -lm4rie -lm4ri -lgivaro -lntl -lgmpxx -lgmp -lm -lstdc++
BUILT_SOURCES =cpucycles.o cpucycles.h

cpucycles.h: cpucycles.o

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; mv cpucycles.o ..; mv cpucycles.h ..; cd ../..)

clean-local:
	-rm cpucycles.h

noinst_PROGRAMS=bench_elimination bench_multiplication

bench_elimination_SOURCES=bench_elimination.cc
bench_elimination_CPPFLAGS=-I.. $(CPPFLAGS)
bench_elimination_LDFLAGS=${MYLDFLAGS} cpucycles.o

bench_multiplication_SOURCES=bench_multiplication.cc
bench_multiplication_CPPFLAGS=-I.. 
bench_multiplication_LDFLAGS=${MYLDFLAGS} cpucycles.o