CPUCYCLES_DIR=./cpucycles-20060326

BENCHCXXFLAGS = $(AM_CFLAGS)  -I.. -I${CPUCYCLES_DIR} @CXXFLAGS@ -DNDEBUG
BENCHLIBADD = -lm4ri cpucycles.o
BENCHLDFLAGS = -Wl,-rpath,../.libs/ ../.libs/libm4rie.so -L${CPUCYCLES_DIR} -Wl,-rpath,${M4RIE_M4RI_PREFIX}/lib

EXTRA_PROGRAMS = bench_multiplication bench_elimination bench_smallops bench_trsm bench_ple

.PHONY: clean dist-clean

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; mv cpucycles.o ..; mv cpucycles.h ..; cd ../..)

cpucycles.h: cpucycles.o

clean-local:
	-rm -f cpucycles.h
	-rm -f $(EXTRA_PROGRAMS)

bench_elimination_SOURCES=bench_elimination.cc cpucycles.h benchmarking.cc
bench_elimination_CXXFLAGS=${BENCHCXXFLAGS}
bench_elimination_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_multiplication_SOURCES=bench_multiplication.cc cpucycles.h benchmarking.cc
bench_multiplication_CXXFLAGS=${BENCHCXXFLAGS}
bench_multiplication_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_smallops_SOURCES=bench_smallops.cc cpucycles.h benchmarking.cc
bench_smallops_CXXFLAGS=${BENCHCXXFLAGS}
bench_smallops_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_trsm_SOURCES=bench_trsm.cc cpucycles.h benchmarking.cc
bench_trsm_CXXFLAGS=${BENCHCXXFLAGS}
bench_trsm_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_ple_SOURCES=bench_ple.cc cpucycles.h benchmarking.cc
bench_ple_CXXFLAGS=${BENCHCXXFLAGS}
bench_ple_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench: cpucycles.o ${EXTRA_PROGRAMS}