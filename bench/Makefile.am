CPUCYCLES_DIR=./cpucycles-20060326

BENCHCXXFLAGS = $(AM_CFLAGS)  -I.. -I${CPUCYCLES_DIR} @CXXFLAGS@ ${M4RIE_GIVARO_CXXFLAGS}
BENCHLIBADD = -lm4ri -lm4rie -lgivaro -lntl -lgmpxx -lgmp -lm -lstdc++ cpucycles.o
BENCHLDFLAGS = ${M4RIE_GIVARO_LDFLAGS} -L${CPUCYCLES_DIR} -Wl,-rpath,${M4RIE_GIVARO_PREFIX}/lib -Wl,-rpath,${M4RIE_M4RI_PREFIX}/lib

EXTRA_PROGRAMS = bench_multiplication bench_elimination bench_smallops

.PHONY: clean dist-clean

cpucycles.h: cpucycles.o

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; mv cpucycles.o ..; mv cpucycles.h ..; cd ../..)

clean-local:
	-rm -f cpucycles.h
	-rm -f $(EXTRA_PROGRAMS)

bench_elimination_SOURCES=bench_elimination.cc cpucycles.h benchmarking.cc
bench_elimination_CXXFLAGS=${BENCHCXXFLAGS}
bench_elimination_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_multiplication_DEPENDS=cpucycles.h cpucycles.o 
bench_multiplication_SOURCES=bench_multiplication.cc cpucycles.h benchmarking.cc
bench_multiplication_CXXFLAGS=${BENCHCXXFLAGS}  -I${CPUCYCLES_DIR}
bench_multiplication_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_smallops_SOURCES=bench_smallops.cc cpucycles.h benchmarking.cc
bench_smallops_CXXFLAGS=${BENCHCXXFLAGS}
bench_smallops_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}
