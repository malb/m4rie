AC_INIT([m4rie],[20200125])

AC_CANONICAL_HOST

AC_CONFIG_SRCDIR(m4rie/newton_john.c)

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE

LT_INIT

AC_PROG_INSTALL

AC_CONFIG_HEADERS(m4rie/config.h)

AC_PROG_CC

AC_SEARCH_LIBS([sqrt], [m], [], [
  AC_MSG_ERROR([unable to find the sqrt() function])
])

if test "$ac_cv_prog_cc_c99" = "no"; then
   AC_MSG_ERROR([C99 support is required but not found.])
fi

# OpenMP support
AC_ARG_ENABLE([openmp],
        AS_HELP_STRING( [--enable-openmp],[add support for OpenMP multicore support.]))

AS_IF([test "x$enable_openmp" = "xyes"], [
   AX_OPENMP()
])
AC_SUBST(OPENMP_CFLAGS)

# Debugging support
AC_ARG_ENABLE(debug, [  --enable-debug          Enable assert() statements for debugging.])

if test "x$enable_debug" = x"yes"; then
   DEBUG_FLAGS="-g"
   AC_SUBST(DEBUG_FLAGS)
else
   AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi

AC_ARG_WITH(m4ri,
        AS_HELP_STRING([--with-m4ri@<:@=VALUE@:>@], [Optional prefix for M4RI library]),[m4rie_config_m4riprefix=$withval])

if test -n $m4rie_config_m4riprefix; then
   M4RIE_M4RI_PREFIX=${m4rie_m4ri_prefix}
   AC_SUBST(M4RIE_M4RI_PREFIX)
   M4RIE_M4RI_LDFLAGS="-L${m4rie_config_m4riprefix}/lib"
   AC_SUBST(M4RIE_M4RI_LDFLAGS)
   M4RIE_M4RI_CFLAGS="-I${m4rie_config_m4riprefix}/include"
   AC_SUBST(M4RIE_M4RI_CFLAGS)
fi

AX_M4RI_CFLAGS()
M4RI_CFLAGS=$ax_cv_m4ri_cflags
AC_SUBST(M4RI_CFLAGS)


RELEASE="AC_PACKAGE_VERSION"
AC_SUBST(RELEASE)

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile tests/Makefile bench/Makefile m4rie.pc])
AC_OUTPUT

