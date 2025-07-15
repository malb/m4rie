AUTOMAKE_OPTIONS = foreign subdir-objects

AM_CFLAGS = ${SIMD_FLAGS} ${DEBUG_FLAGS} ${M4RI_CFLAGS} -I${top_srcdir}
LDADD = ${top_builddir}/libm4rie.la @M4RI_LIBS@
AM_LDFLAGS = ${M4RI_LIBS} -no-install

EXTRA_DIST = testing.h

TESTS = test_trsm test_elimination test_multiplication test_smallops test_ple test_transpose
check_PROGRAMS = ${TESTS}

all: ${TESTS}
