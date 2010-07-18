AUTOMAKE_OPTIONS = gnu

AM_CFLAGS=${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm4rie.la

libm4rie_la_SOURCES = src/finite_field.c src/gf2e_matrix.c src/travolta.c

pkgincludesubdir = $(includedir)/m4rie
pkgincludesub_HEADERS = src/finite_field.h  src/gf2e_matrix.h  src/m4rie.h  src/m4ri_functions.h  src/travolta.h

libm4rie_la_LDFLAGS = -release 0.0.19780616

check_PROGRAMS=test_elimination
test_elimination_SOURCES=tests/test_elimination.cc
test_elimination_LDFLAGS= -lm4rie -lm4ri -lgivaro -lntl -lgmpxx -lgmp -lm -lstdc++

TESTS = test_elimination
