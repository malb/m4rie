AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS=${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS} ${M4RIE_M4RI_CFLAGS} ${M4RI_CFLAGS}
EXTRA_DIST = m4rie/Doxyfile

lib_LTLIBRARIES = libm4rie.la

libm4rie_la_SOURCES = m4rie/gf2e.c \
	m4rie/mzed.c \
	m4rie/newton_john.c \
	m4rie/echelonform.c \
	m4rie/strassen.c \
	m4rie/mzd_slice.c \
	m4rie/mzd_poly.c \
	m4rie/mzd_ptr.c \
	m4rie/karatsuba.c \
	m4rie/blm.c \
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
	m4rie/blm.h \
	m4rie/trsm.h \
	m4rie/ple.h \
	m4rie/permutation.h \
	m4rie/conversion.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = m4rie.pc

libm4rie_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined ${M4RIE_M4RI_LDFLAGS}
libm4rie_la_LIBADD = -lm4ri

SUBDIRS = . tests
DIST_SUBDIRS = $(SUBDIRS) bench

