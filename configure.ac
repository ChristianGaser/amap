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
      
dnl Build only static libs by default

AC_FUNC_FORK
AC_CHECK_HEADERS(values.h malloc.h)
AC_CHECK_FUNCS(finite isfinite mkstemp)

AC_CHECK_LIB([m], [sqrt], LDFLAGS="$LDFLAGS -lm")

LDFLAGS="$LDFLAGS"

AC_CHECK_HEADERS()

AC_CONFIG_FILES([Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_OUTPUT
