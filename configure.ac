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

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile bench/Makefile])
AC_OUTPUT

