#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pdp3, CURRENT)

dnl m4_pattern_allow([^AC_])

dnl AC_CONFIG_AUX_DIR([acdir])
dnl AC_CONFIG_MACRO_DIR([acdir])

# AC_CANONICAL_SYSTEM

AC_SUBST(EXPERIMENTAL_OPTION)
AC_SUBST(SPEEDUP)
AC_SUBST(DEBUG)
AC_SUBST(SINGLETHREAD)
AC_SUBST(HDF5)
AC_SUBST(CFLAGS_DEFAULT)
AC_SUBST(CFLAGS_ADDITIONAL)
AC_SUBST(DOXYGEN_FORMATS)
AC_SUBST(HDF5_OPTION)
AC_SUBST(DEBUG_OPTION)

AC_SUBST(MKDIR)
MKDIR=mkdir

AC_SUBST(RESULT_DIR)

AC_SUBST(RELEASE)
RELEASE=CURRENT

AC_ARG_ENABLE([experimental], [AC_HELP_STRING([--enable-experimental], [Enable experimental features (only for debug/development)])], [ENABLE_EXPERIMENTAL="$enableval"], [ENABLE_EXPERIMENTAL=no])
AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug], [Enable debug (optimized for GDB)])], [DEBUG="$enableval"], [DEBUG=no])
AC_ARG_ENABLE([speedup], [AC_HELP_STRING([--enable-speedup], [Enable speedup])], [SPEEDUP="$enableval"], [SPEEDUP=no])
AC_ARG_ENABLE([singlethread], [AC_HELP_STRING([--enable-singlethread], [Build application to run in single-thread mode (used in debug by default)])], [SINGLETHREAD="$enableval"], [SINGLETHREAD=${DEBUG}])

AC_ARG_ENABLE([hdf5], [AC_HELP_STRING([--enable-hdf5], [Build application with hdf5_serial library support])], [HDF5="$enableval"], [HDF5=no])

AC_ARG_WITH([doc-formats], [AC_HELP_STRING([--with-doc-formats], [Set list of documentation formats to be generated])], [DOXYGEN_FORMATS="$withval"], [DOXYGEN_FORMATS="latex html rtf"])
AC_ARG_WITH([result-dir], [AC_HELP_STRING([--with-result-dir], [Set default result directory])], [RESULT_DIR="$withval"], [RESULT_DIR="pdp3_result"])

AC_PROG_CXXCPP
AC_LANG(C++)

dnl # Checks for libraries.
if test x$SINGLETHREAD = xno; then
   AC_OPENMP
fi

if test x$HDF5 = xyes; then
   AC_HAVE_LIBRARY(hdf5_serial, , AC_MSG_ERROR([No such library hdf5_serial.]))
fi

dnl # Checks for header files.

dnl # Checks for typedefs, structures, and compiler characteristics.

# set default and customized options
CFLAGS_DEFAULT='-m64 -mcmodel=medium -std=c++14'

# experimental option
if test x$ENABLE_EXPERIMENTAL = xyes; then
   EXPERIMENTAL_OPTION="-DEXPERIMENTAL"
fi

# speedup
if test x$SPEEDUP = xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -ffast-math -O3 -funroll-loops -Wno-write-strings -freciprocal-math -fcx-fortran-rules"
else
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -O2"
fi

# singlethread
if test x$SINGLETHREAD = xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -Wno-unknown-pragmas"
else
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -fopenmp"
fi

# debug
if test x$DEBUG = xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -O0 -Wall -Wextra -ggdb3 -fvar-tracking -ggnu-pubnames -pedantic -time -ftree-vectorizer-verbose=7 -DDEBUG -DPERF_DEBUG"
  DEBUG_OPTION=true
else
  DEBUG_OPTION=false
fi

if test x$HDF5 = xyes; then
   HDF5_OPTION=true
else
   HDF5_OPTION=false
fi

# Checks for library functions.
AC_OUTPUT(Makefile)
AC_OUTPUT(parameters.xml)
