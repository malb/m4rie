AUTOMAKE_OPTIONS = gnu

AM_CFLAGS=${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS} ${M4RIE_M4RI_CFLAGS}

lib_LTLIBRARIES = libm4rie.la

libm4rie_la_SOURCES = src/finite_field.c src/gf2e_matrix.c src/travolta.c src/echelonform.c src/strassen.c src/bitslice.c

pkgincludesubdir = $(includedir)/m4rie
pkgincludesub_HEADERS = src/finite_field.h  src/gf2e_matrix.h  src/m4rie.h  src/m4ri_functions.h  src/travolta.h src/echelonform.h src/strassen.h src/bitslice.h gf2e_cxx/finite_field_givaro.h 

libm4rie_la_LDFLAGS = -release 0.0.20110801 -no-undefined ${M4RIE_M4RI_LDFLAGS}
libm4rie_la_LIBADD = -lm4ri

#testing
TESTCXXFLAGS = $(AM_CFLAGS) @CXXFLAGS@ ${M4RIE_M4RI_CFLAGS} ${M4RIE_GIVARO_CXXFLAGS}
TESTLIBADD = -lm4ri -lm4rie -lgivaro -lntl -lgmpxx -lgmp -lm -lstdc++
TESTLDFLAGS = ${M4RIE_GIVARO_LDFLAGS}

check_PROGRAMS=test_elimination test_multiplication test_smallops
test_elimination_SOURCES=tests/test_elimination.cc
test_elimination_LDADD= ${TESTLDFLAGS} ${TESTLIBADD}
test_elimination_CXXFLAGS=${TESTCXXFLAGS}

test_multiplication_SOURCES=tests/test_multiplication.cc
test_multiplication_LDADD=${TESTLDFLAGS} ${TESTLIBADD}
test_multiplication_CXXFLAGS=${TESTCXXFLAGS}

test_smallops_SOURCES=tests/test_smallops.cc
test_smallops_LDADD=${TESTLDFLAGS} ${TESTLIBADD}
test_smallops_CXXFLAGS=${TESTCXXFLAGS}

TESTS = test_elimination test_multiplication test_smallops

# benchmarketing

SUBDIRS = . bench

clean-local:
	(cd tests; make clean; cd ..)