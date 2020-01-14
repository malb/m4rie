AUTOMAKE_OPTIONS = foreign subdir-objects

AM_CFLAGS = ${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS} ${M4RIE_M4RI_CFLAGS} ${M4RI_CFLAGS} -I${top_srcdir} @CFLAGS@
LDADD = ${top_builddir}/libm4rie.la -lm4ri -lm
LDFLAGS = ${M4RIE_M4RI_LDFLAGS} -no-install

EXTRA_DIST = testing.h

TESTS = test_trsm test_elimination test_multiplication test_smallops test_ple
check_PROGRAMS = ${TESTS}

all: ${TESTS}
