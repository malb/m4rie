AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS=${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS} ${M4RIE_M4RI_CFLAGS} ${M4RI_CFLAGS}

lib_LTLIBRARIES = libm4rie.la

libm4rie_la_SOURCES = src/gf2e.c \
	src/mzed.c \
	src/newton_john.c \
	src/echelonform.c \
	src/strassen.c \
	src/mzd_slice.c \
	src/mzd_poly.c \
	src/trsm.c \
	src/ple.c \
	src/conversion.c

pkgincludesubdir = $(includedir)/m4rie

pkgincludesub_HEADERS = src/gf2e.h \
	src/mzed.h \
	src/m4rie.h \
	src/m4ri_functions.h \
	src/newton_john.h \
	src/echelonform.h \
	src/strassen.h \
	src/mzd_slice.h \
	src/mzd_poly.h \
	src/trsm.h \
	src/ple.h \
	src/permutation.h \
	src/conversion.h \
	gf2e_cxx/finite_field_givaro.h

libm4rie_la_LDFLAGS = -release 0.0.20120415 -no-undefined ${M4RIE_M4RI_LDFLAGS}
libm4rie_la_LIBADD = -lm4ri

#testing
TESTCXXFLAGS = ${AM_CFLAGS} @CXXFLAGS@ ${M4RIE_M4RI_CFLAGS} ${M4RIE_GIVARO_CXXFLAGS} ${M4RI_CFLAGS}
TESTLIBADD = -lm4ri -lm4rie -lgivaro -lgmpxx -lgmp -lm -lstdc++
TESTLDADD = ${M4RIE_GIVARO_LDFLAGS}
TESTLDFLAGS =  -Wl,-rpath,.libs/

check_PROGRAMS=test_trsm test_elimination test_multiplication test_smallops test_ple
test_elimination_SOURCES=tests/test_elimination.cc
test_elimination_LDADD= ${TESTLDADD} ${TESTLIBADD}
test_elimination_LDFLAGS= ${TESTLDFLAGS}
test_elimination_CXXFLAGS=${TESTCXXFLAGS}

test_multiplication_SOURCES=tests/test_multiplication.cc
test_multiplication_LDADD=${TESTLDADD} ${TESTLIBADD}
test_multiplication_LDFLAGS= ${TESTLDFLAGS}
test_multiplication_CXXFLAGS=${TESTCXXFLAGS}

test_smallops_SOURCES=tests/test_smallops.cc
test_smallops_LDADD=${TESTLDADD} ${TESTLIBADD}
test_smallops_LDFLAGS= ${TESTLDFLAGS}
test_smallops_CXXFLAGS=${TESTCXXFLAGS}

test_trsm_SOURCES=tests/test_trsm.cc
test_trsm_LDADD=${TESTLDADD} ${TESTLIBADD}
test_trsm_LDFLAGS= ${TESTLDFLAGS}
test_trsm_CXXFLAGS=${TESTCXXFLAGS}

test_ple_SOURCES=tests/test_ple.cc
test_ple_LDADD=${TESTLDADD} ${TESTLIBADD}
test_ple_LDFLAGS= ${TESTLDFLAGS}
test_ple_CXXFLAGS=${TESTCXXFLAGS}


TESTS=test_trsm test_elimination test_multiplication test_smallops test_ple

# benchmarketing

SUBDIRS = . bench

clean-local:
	(cd tests; make clean; cd ..)
