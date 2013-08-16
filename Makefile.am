AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS=${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS} ${M4RIE_M4RI_CFLAGS} ${M4RI_CFLAGS}

lib_LTLIBRARIES = libm4rie.la

libm4rie_la_SOURCES = m4rie/gf2e.c \
	m4rie/mzed.c \
	m4rie/newton_john.c \
	m4rie/echelonform.c \
	m4rie/strassen.c \
	m4rie/mzd_slice.c \
	m4rie/mzd_poly.c \
	m4rie/mzd_ptr.c \
	m4rie/trsm.c \
	m4rie/ple.c \
	m4rie/conversion.c \
	m4rie/conversion_slice8.c \
	m4rie/conversion_slice16.c \
	m4rie/conversion_cling8.c \
	m4rie/conversion_cling16.c \
	m4rie/mzd_slice_intro.inl \
	m4rie/mzd_slice_outro.inl \
	m4rie/mzed_intro.inl \
	m4rie/mzed_outro.inl \
	m4rie/trsm.inl

pkgincludesubdir = $(includedir)/m4rie

pkgincludesub_HEADERS = m4rie/gf2x.h \
	m4rie/gf2e.h \
	m4rie/mzed.h \
	m4rie/m4rie.h \
	m4rie/m4ri_functions.h \
	m4rie/newton_john.h \
	m4rie/echelonform.h \
	m4rie/strassen.h \
	m4rie/mzd_slice.h \
	m4rie/mzd_poly.h \
	m4rie/mzd_ptr.h \
	m4rie/trsm.h \
	m4rie/ple.h \
	m4rie/permutation.h \
	m4rie/conversion.h \
	gf2e_cxx/finite_field_givaro.h

libm4rie_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined ${M4RIE_M4RI_LDFLAGS}
libm4rie_la_LIBADD = -lm4ri

#testing
TESTCXXFLAGS = ${AM_CFLAGS} @CXXFLAGS@ ${M4RIE_M4RI_CFLAGS} ${M4RI_CFLAGS} -I./tests
TESTLIBADD = -lm4ri -lm4rie 
TESTLDADD = 
TESTLDFLAGS = ${M4RIE_M4RI_LDFLAGS}

check_PROGRAMS = test_trsm test_elimination test_multiplication test_smallops test_ple
test_elimination_SOURCES = tests/test_elimination.cc tests/testing.h
test_elimination_LDADD =  ${TESTLDADD} ${TESTLIBADD}
test_elimination_LDFLAGS = ${TESTLDFLAGS}
test_elimination_CXXFLAGS = ${TESTCXXFLAGS}

test_multiplication_SOURCES = tests/test_multiplication.cc tests/testing.h
test_multiplication_LDADD = ${TESTLDADD} ${TESTLIBADD}
test_multiplication_LDFLAGS = ${TESTLDFLAGS}
test_multiplication_CXXFLAGS = ${TESTCXXFLAGS}

test_smallops_SOURCES = tests/test_smallops.cc tests/testing.h
test_smallops_LDADD = ${TESTLDADD} ${TESTLIBADD}
test_smallops_LDFLAGS = ${TESTLDFLAGS}
test_smallops_CXXFLAGS = ${TESTCXXFLAGS}

test_trsm_SOURCES = tests/test_trsm.cc tests/testing.h
test_trsm_LDADD = ${TESTLDADD} ${TESTLIBADD}
test_trsm_LDFLAGS = ${TESTLDFLAGS}
test_trsm_CXXFLAGS = ${TESTCXXFLAGS}

test_ple_SOURCES = tests/test_ple.cc tests/testing.h
test_ple_LDADD = ${TESTLDADD} ${TESTLIBADD}
test_ple_LDFLAGS = ${TESTLDFLAGS}
test_ple_CXXFLAGS = ${TESTCXXFLAGS}

TESTS = test_trsm test_elimination test_multiplication test_smallops test_ple

# benchmarketing

SUBDIRS = . bench

clean-local:
	(cd tests; make clean; cd ..)
