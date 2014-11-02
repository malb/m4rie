CPUCYCLES_DIR=./cpucycles-20060326

CFLAGS = $(AM_CFLAGS)  -I.. -I${CPUCYCLES_DIR} @CXXFLAGS@ -DNDEBUG ${M4RIE_M4RI_CFLAGS}
LDADD = ${top_builddir}/libm4rie.la -lm4ri -lm cpucycles.o benchmarking.o
LDFLAGS = -Wl,-rpath,../.libs/ ../.libs/libm4rie.so -L${CPUCYCLES_DIR} ${M4RIE_M4RI_LDFLAGS} -Wl,-rpath,${M4RIE_M4RI_PREFIX}/lib

EXTRA_PROGRAMS = bench_multiplication bench_elimination bench_smallops bench_trsm bench_ple

.PHONY: clean dist-clean

cpucycles.o:
	(cd $(CPUCYCLES_DIR); sh do; mv cpucycles.o ..; mv cpucycles.h ..; cd ../..)

cpucycles.h: cpucycles.o

clean-local:
	-rm -f cpucycles.h
	-rm -f $(EXTRA_PROGRAMS)

all: ${EXTRA_PROGRAMS}
