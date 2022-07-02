dnl $Id: configure.in 214 2019-01-31 11:20:14Z gaser $
dnl Process this file with autoconf to produce a configure script.

AC_INIT(amap,0.1,Christian Gaser <christian.gaser@uni-jena.de>)
AM_INIT_AUTOMAKE

AC_REVISION($Revision: 214 $)

AM_CONFIG_HEADER([config.h])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
dnl AC_PROG_LN_S
dnl AC_PROG_MAKE_SET
dnl AC_PROG_RANLIB
LT_INIT

AC_ENABLE_STATIC
AC_DISABLE_SHARED

smr_WITH_BUILD_PATH

AC_OPENMP
AC_CHECK_HEADERS([omp.h],[
  AC_DEFINE(OPENMP,1,[Define to enable OpenMP.])
    CFLAGS="$CFLAGS $OPENMP_CFLAGS"
    CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
  CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
  LDFLAGS="$LDFLAGS $OPENMP_CFLAGS"
])

AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp],[enable OpenMP support [default=yes]])], enable_openmp=$enableval, enable_openmp=yes)
if test "$enable_openmp" = "yes"; then
  AC_DEFINE(OPENMP,1,[Define to enable OpenMP.])
else
  AC_DEFINE(OPENMP,0,[Define to enable OpenMP.])
fi
AM_CONDITIONAL([OPENMP], [test "$enable_openmp" = yes])


dnl Build only static libs by default

AC_FUNC_FORK
AC_CHECK_HEADERS(values.h malloc.h)
AC_CHECK_FUNCS(finite isfinite mkstemp)

AC_CHECK_LIB([m], [sqrt], LDFLAGS="$LDFLAGS -lm")

AC_CHECK_HEADERS()

AC_CONFIG_FILES([Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_OUTPUT
