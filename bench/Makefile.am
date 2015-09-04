# TOPSRCDIR allows to compile this testsuite for a clone (and different revision) of this
# repository. If the environment variable topsrcdir is the top source directory of the clone,
# compile this testsuite as:
#
# make TOPSRCDIR="$topsrcdir"
#
# Finally, if you also leave out TOPSRCDIR and just run 'make' then the testsuite
# is compiled against the current source tree.

AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = $(srcdir)/..
TOPBUILDDIR = $(builddir)/..
DEFINES =
# include TOPBUILDIR for m4ri_config.h
AM_CFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR) -D_XOPEN_SOURCE=600 @CFLAGS@ $(DEFINES) @OPENMP_CFLAGS@ ${M4RIE_M4RI_CFLAGS}

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) cpucycles.o -lm4ri -lm4rie -lm ${M4RIE_M4RI_LDFLAGS} -Wl,-rpath,${M4RIE_M4RI_PREFIX}/lib

BENCH = bench_multiplication \
        bench_elimination \
        bench_smallops \
        bench_trsm \
        bench_ple


CPUCYCLES_DIR = cpucycles-20060326

CPUCYCLES_EXTRA_DIST = $(srcdir)/cpucycles-20060326/*.c \
                       $(srcdir)/cpucycles-20060326/*.h \
                       $(srcdir)/cpucycles-20060326/compile \
                       $(srcdir)/cpucycles-20060326/do

EXTRA_DIST = $(CPUCYCLES_EXTRA_DIST)

bin_PROGRAMS = $(BENCH)

bench_elimination_SOURCES = bench_elimination.c benchmarking.c benchmarking.h
bench_multiplication_SOURCES = bench_multiplication.c benchmarking.c benchmarking.h
bench_ple_SOURCES = bench_ple.c benchmarking.c benchmarking.h
bench_smallops_SOURCES = bench_smallops.c benchmarking.c benchmarking.h
bench_trsm_SOURCES = bench_trsm.c benchmarking.c benchmarking.h

BUILT_SOURCES = cpucycles.h

cpucycles.h: cpucycles.o

cpucycles.o:
	 (if [ $(srcdir) != $(builddir) ]; then \
		cp -r $(srcdir)/$(CPUCYCLES_DIR) $(builddir); \
	        chmod +w -R $(CPUCYCLES_DIR)/; \
         fi; \
         cd $(CPUCYCLES_DIR); \
         sh do; \
         cp cpucycles.o ..; \
         cp cpucycles.h ..; \
         cd ..; \
	 if [ $(srcdir) != $(builddir) ]; then \
		rm -r $(CPUCYCLES_DIR); \
         fi)

distclean-local:
	-rm -f cpucycles.h
