#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
# AC_INIT(pdp3, CURRENT)

AC_INIT(pdp3, esyscmd(git log --pretty=%h -1 | tr -d \\n), cosmonaut.ok@zoho.com)

AC_CONFIG_HEADERS([include/config.h])
AC_OUTPUT

m4_pattern_allow([^AC_])

dnl AC_CONFIG_AUX_DIR([acdir])
dnl AC_CONFIG_MACRO_DIR([acdir])

# AC_CANONICAL_SYSTEM

AC_SUBST(EXPERIMENTAL_OPTION)
AC_SUBST(LEGACY_OPTION)
AC_SUBST(IEEE)
AC_SUBST(DEBUG)
AC_SUBST(TESTMODE)
AC_SUBST(SINGLETHREAD)
AC_SUBST(HDF5)
AC_SUBST(CFLAGS_DEFAULT)
AC_SUBST(CFLAGS_ADDITIONAL)
AC_SUBST(LDFLAGS_ADDITIONAL)
AC_SUBST(DOXYGEN_FORMATS)
AC_SUBST(HDF5_OPTION)
AC_SUBST(PUSHER)
AC_SUBST(OMP_DYNAMIC)
AC_SUBST(DEBUG_OPTION)
AC_SUBST(PROFILER)
AC_SUBST(HDF5_URL)
AC_SUBST(HDF5_VERSION)

AC_SUBST(MKDIR)
MKDIR=mkdir

AC_SUBST(RELEASE)
RELEASE=${PACKAGE_VERSION}

HDF5_URL="https://www.hdfgroup.org/package/source-bzip/?wpdmdl=12594&refresh=5bbc7778635b21539078008"
HDF5_VERSION="hdf5-1.10.3"

AC_ARG_ENABLE([experimental], [AC_HELP_STRING([--enable-experimental], [Enable experimental features (only for debug/development)])], [ENABLE_EXPERIMENTAL="$enableval"], [ENABLE_EXPERIMENTAL=no])
AC_ARG_ENABLE([legacy], [AC_HELP_STRING([--enable-legacy], [Use legacy instead of current features])], [ENABLE_LEGACY="$enableval"], [ENABLE_LEGACY=no])
AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug], [Enable debug (optimized for GDB)])], [DEBUG="$enableval"], [DEBUG=no])
AC_ARG_ENABLE([testmode], [AC_HELP_STRING([--enable-testmode], [Enable test mode with reproducible random numbers to compare result (USE ONLY FOR FUNCTIONAL OR OTHER TESTS)])], [TESTMODE="$enableval"], [TESTMODE=no])
AC_ARG_ENABLE([profiler], [AC_HELP_STRING([--enable-profiler], [Enable profiler support for performance debug])], [PROFILER="$enableval"], [PROFILER=no])
AC_ARG_ENABLE([ieee], [AC_HELP_STRING([--enable-ieee], [Keep IEEE compliant rounding and disable hardware acceleration. Decrease calculation speed])], [IEEE="$enableval"], [IEEE=no])
AC_ARG_ENABLE([singlethread], [AC_HELP_STRING([--enable-singlethread], [Build application to run in single-thread mode (used in debug by default)])], [SINGLETHREAD="$enableval"], [SINGLETHREAD=${DEBUG}])
AC_ARG_ENABLE([omp_dynamic], [AC_HELP_STRING([--enable-omp-dynamic], [Enable openmp dynamic thread creation])], [OMP_DYNAMIC="$enableval"], [OMP_DYNAMIC=no])

AC_ARG_ENABLE([hdf5], [AC_HELP_STRING([--enable-hdf5], [Build application with hdf5 library support])], [HDF5="$enableval"], [HDF5=no])

AC_ARG_WITH([doc-formats], [AC_HELP_STRING([--with-doc-formats], [Set list of documentation formats to be generated])], [DOXYGEN_FORMATS="$withval"], [DOXYGEN_FORMATS="latex html rtf"])
AC_ARG_WITH([pusher], [AC_HELP_STRING([--with-pusher], [Set particles pusher (supported: boris-adaptive (default), boris, boris-relativistic, vay, hc (n/a), rrll (n/a))])], [PUSHER="$withval"], [PUSHER="boris-adaptive"])



AC_ARG_ENABLE([temperature-weighted], [AC_HELP_STRING([--enable-temperature-weighted], [Enable temperature calculation algorithm, based on weighting, instead of classic calculations])], [TEMPERATURE_WEIGHTED="$enableval"], [TEMPERATURE_WEIGHTED=no])

AC_ARG_WITH([temperature-postproc], [AC_HELP_STRING([--with-temperature-postproc], [Set temperature postprocessing to: none or bilinear])], [TEMPERATURE_POSTPROC="$withval"], [TEMPERATURE_POSTPROC="none"])

AC_PROG_CXXCPP
AC_LANG(C++)

dnl # Checks for libraries.
if test x$SINGLETHREAD = xno; then
   AC_OPENMP
fi

if test x$HDF5 = xyes; then
   AC_HAVE_LIBRARY(hdf5, , AC_MSG_WARN([No such library hdf5. You can build it locally after configuration with `make hdf5-install` before build and it includes to build paths automatically]))
fi

dnl # Checks for header files.

dnl # Checks for typedefs, structures, and compiler characteristics.

# set default and customized options
CFLAGS_DEFAULT="-m64 -mcmodel=medium -std=c++14"

# experimental option
if test x$ENABLE_EXPERIMENTAL = xyes; then
   AC_DEFINE([EXPERIMENTAL], [true], [Use some experimental features])
fi

# legacy option
if test x$ENABLE_LEGACY = xyes; then
   AC_DEFINE([LEGACY], [true], [Use legacy features instead of current recommended])
fi

# speedup
if test x$IEEE != xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -ffast-math -O3 -march=native -funroll-loops -Wno-write-strings -freciprocal-math -fcx-fortran-rules -msse -mavx"
  AC_DEFINE([SPEEDUP], [true], [Directive to increase calculation speed])
else
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -O2"
fi

# singlethread
if test x$SINGLETHREAD = xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -Wno-unknown-pragmas"
  AC_DEFINE([SINGLETHREAD], [true], [Debug option])
else
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -fopenmp"
fi

# debug
if test x$DEBUG = xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -O0 -Wall -Wextra -ggdb3 -fvar-tracking -ggnu-pubnames -pedantic -time -ftree-vectorizer-verbose=7"
  DEBUG_OPTION=true
  AC_DEFINE([DEBUG], [true], [Use debug option])
else
  DEBUG_OPTION=false
fi

# debug
if test x$TESTMODE = xyes; then
  AC_DEFINE([TESTMODE], [true], [Use testmode. Only for functional tests with reproducible random numbers])
else
  DEBUG_OPTION=false
fi

# profiler
if test x$PROFILER = xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -pg -no-pie"
  LDFLAGS_ADDITIONAL="$LDFLAGS_ADDITIONAL -pg -no-pie"
  AC_DEFINE([PROFILER], [true], [Performance debug option])
fi

if test x$HDF5 = xyes; then
   HDF5_OPTION=true
   AC_DEFINE([USE_HDF5], [true], [Use HDF5])
else
   HDF5_OPTION=false
fi

# define a pusher
if test x$PUSHER == xboris-adaptive; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_ADAPTIVE], [true], [Performance debug option])
   AC_DEFINE([REL_LIMIT], [5e7], [Relativistit calculations limit (for boris-adaptive pusher)])
elif test x$PUSHER == xboris; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_CLASSIC], [true], [Performance debug option])
elif test x$PUSHER == xboris-relativistic; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_RELATIVISTIC], [true], [Performance debug option])
elif test x$PUSHER == xvay; then
   AC_DEFINE_UNQUOTED([PUSHER_VAY], [true], [Performance debug option])
else
   AC_MSG_ERROR([Pusher is not supported.])
fi

# define a pusher
if test x$TEMPERATURE_POSTPROC == xnone; then
   AC_DEFINE_UNQUOTED([TEMPERATURE_POSTPROC_NONE], [true], [temperature postprocessing])
elif test x$TEMPERATURE_POSTPROC == xbilinear; then
   AC_DEFINE_UNQUOTED([TEMPERATURE_POSTPROC_BILINEAR], [true], [temperature postprocessing])
dnl elif test x$TEMPERATURE_POSTPROC == xbicubic; then # FIXME: bicubic still not working correctly
dnl    AC_DEFINE_UNQUOTED([TEMPERATURE_POSTPROC_BICUBIC], [true], [temperature postprocessing])
else
   AC_MSG_ERROR([Temperature postprocessing method is not supported.])
fi

if test x$OMP_DYNAMIC == xyes; then
   AC_DEFINE([OPENMP_DYNAMIC_THREADS], [true], [Set OpenMP threads to be created dynamically])
fi

if test x$TEMPERATURE_WEIGHTED == xyes; then
   AC_DEFINE([TEMPERATURE_WEIGHTED], [true], [Enable new temperature calculation algorithm])
fi

# Checks for library functions.
AC_OUTPUT(Makefile)
AC_OUTPUT(parameters.xml)
AC_OUTPUT(test/functional/parameters.xml.tmpl)
AC_OUTPUT(test/functional/parameters_ext.xml.tmpl)
AC_OUTPUT(test/functional/parameters_regression.xml.tmpl)
