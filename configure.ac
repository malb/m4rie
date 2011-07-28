AC_INIT(libm4rie,1.0.1)

dnl AC_PREREQ([2.64])

AC_CANONICAL_HOST

AC_CONFIG_SRCDIR(src/travolta.c)

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE

AC_PROG_LIBTOOL

AC_PROG_INSTALL

AC_CONFIG_HEADERS(src/config.h)

AC_PROG_CXX

m4_include([m4/ax_check_compiler_flags.m4])dnl
m4_include([m4/ax_openmp.m4])dnl

AC_PROG_CC_C99()


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

AC_ARG_WITH(givaro,
        AS_HELP_STRING([--with-givaro@<:@=VALUE@:>@], [Optional prefix for Givaro library (for testing).]),[m4rie_config_givaroprefix=$withval])
if test -n $m4rie_config_givaroprefix; then
   M4RIE_GIVARO_PREFIX=${m4rie_givaro_prefix}
   AC_SUBST(M4RIE_GIVARO_PREFIX)
   M4RIE_GIVARO_LDFLAGS="-L${m4rie_config_givaroprefix}/lib"
   AC_SUBST(M4RIE_GIVARO_LDFLAGS)
   M4RIE_GIVARO_CXXFLAGS="-I${m4rie_config_givaroprefix}/include"
   AC_SUBST(M4RIE_GIVARO_CXXFLAGS)
fi

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile bench/Makefile])
AC_OUTPUT

