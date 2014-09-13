CPUCYCLES_DIR=./cpucycles-20060326

BENCHCFLAGS = $(AM_CFLAGS)  -I.. -I${CPUCYCLES_DIR} @CXXFLAGS@ -DNDEBUG ${M4RIE_M4RI_CFLAGS}
BENCHLIBADD = -lm4ri -lm cpucycles.o
BENCHLDFLAGS = -Wl,-rpath,../.libs/ ../.libs/libm4rie.so -L${CPUCYCLES_DIR} ${M4RIE_M4RI_LDFLAGS} -Wl,-rpath,${M4RIE_M4RI_PREFIX}/lib

EXTRA_PROGRAMS = bench_multiplication bench_elimination bench_smallops bench_trsm bench_ple

.PHONY: clean dist-clean

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; mv cpucycles.o ..; mv cpucycles.h ..; cd ../..)

cpucycles.h: cpucycles.o

clean-local:
	-rm -f cpucycles.h
	-rm -f $(EXTRA_PROGRAMS)

bench_elimination_SOURCES=bench_elimination.c cpucycles.h benchmarking.c
bench_elimination_CFLAGS=${BENCHCFLAGS}
bench_elimination_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_multiplication_SOURCES=bench_multiplication.c cpucycles.h benchmarking.c
bench_multiplication_CFLAGS=${BENCHCFLAGS}
bench_multiplication_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_smallops_SOURCES=bench_smallops.c cpucycles.h benchmarking.c
bench_smallops_CFLAGS=${BENCHCFLAGS}
bench_smallops_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_trsm_SOURCES=bench_trsm.c cpucycles.h benchmarking.c
bench_trsm_CFLAGS=${BENCHCFLAGS}
bench_trsm_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench_ple_SOURCES=bench_ple.c cpucycles.h benchmarking.c
bench_ple_CFLAGS=${BENCHCFLAGS}
bench_ple_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

bench: cpucycles.o ${EXTRA_PROGRAMS}
