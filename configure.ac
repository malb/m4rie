AC_INIT([m4rie],[20250128])

# See http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
LT_VERSION=1:0:0
AC_SUBST(LT_VERSION)

AC_CANONICAL_HOST

AC_CONFIG_SRCDIR(m4rie/newton_john.c)

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE

LT_INIT

AC_PROG_INSTALL

AC_CONFIG_HEADERS(m4rie/config.h)

AC_PROG_CC

# libm is only used directly in bench/benchmarking.c
AC_SEARCH_LIBS([sqrt], [m], [], [
  AC_MSG_ERROR([unable to find the sqrt() function])
])

if test "$ac_cv_prog_cc_c99" = "no"; then
   AC_MSG_ERROR([C99 support is required but not found.])
fi

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

AS_IF([test -n "$m4rie_config_m4riprefix"], [
   M4RIE_M4RI_PREFIX=${m4rie_m4ri_prefix}
   AC_SUBST(M4RIE_M4RI_PREFIX)

   M4RI_CFLAGS="-I${m4rie_config_m4riprefix}/include"
   M4RI_LIBS="-L${m4rie_config_m4riprefix}/lib -lm4ri -lm"
], [
  # Don't use pkgconfig if the user passed --with-m4ri=$path,
  # because in that case, we want to use the copy from $path
  # and not the one that pkgconfig knows about.
  PKG_CHECK_MODULES([M4RI], [m4ri >= 20240729], [], [
    # pkgconfig failed, but we can try the old-fashioned way
    AC_SEARCH_LIBS([mzd_info],[m4ri],[
      # Current versions of m4ri.pc have "Libs: ... -lm". It's not
      # 100% clear if libm is supposed to be part of its interface,
      # but m4ri.h includes math.h, so we better play it safe.
      M4RI_LIBS="-lm4ri -lm"
    ], [
      AC_MSG_ERROR([required library m4ri not found])
    ])
  ])
])

AC_SUBST(M4RI_CFLAGS)
AC_SUBST(M4RI_LIBS)


RELEASE="AC_PACKAGE_VERSION"
AC_SUBST(RELEASE)

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile tests/Makefile bench/Makefile m4rie.pc])
AC_OUTPUT

