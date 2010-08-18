CPUCYCLES_DIR=./cpucycles-20060326
MYLDFLAGS=-L../.libs -lm4rie -lm4ri -lgivaro -lntl -lgmpxx -lgmp -lm -lstdc++

cpucycles.h: cpucycles.o

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; mv cpucycles.o ..; mv cpucycles.h ..; cd ../..)

EXTRA_PROGRAMS = bench_multiplication bench_elimination

clean-local:
	-rm -f cpucycles.h
	-rm -f $(EXTRA_PROGRAMS)

bench_elimination_SOURCES=bench_elimination.cc
bench_elimination_CPPFLAGS=-I.. $(CPPFLAGS)
bench_elimination_LDFLAGS=${MYLDFLAGS} cpucycles.o

bench_multiplication_SOURCES=bench_multiplication.cc
bench_multiplication_CPPFLAGS=-I.. $(CPPFLAGS)
bench_multiplication_LDFLAGS=${MYLDFLAGS} cpucycles.o

