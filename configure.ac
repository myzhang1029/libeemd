#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([libeemd], [1.4.1], [https://bitbucket.org/luukko/libeemd/issues])
AC_CONFIG_SRCDIR([src/eemd.c])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
LT_INIT([win32-dll])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S

# Optionally build with OpenMP support
AC_ARG_ENABLE([openmp],
    AS_HELP_STRING([--disable-openmp], [Do not build with OpenMP parallelism.]))

# Checks for libraries.
AC_CHECK_LIB([gsl], [gsl_strerror], [], [
              AC_MSG_ERROR([Cannot find libgsl. Try setting LDFLAGS and CFLAGS.])
])

# Checks for header files.
AC_CHECK_HEADERS([
        gsl/gsl_statistics_double.h gsl/gsl_rng.h gsl/gsl_randist.h
        gsl/gsl_vector.h gsl/gsl_linalg.h gsl/gsl_poly.h
    ], [], [
        AC_MSG_ERROR([Cannot find gsl headers. Try setting CFLAGS.])
])

# Enable OpenMP if found
AC_OPENMP
AS_IF([test "x${ac_cv_prog_c_openmp}" == "xunsupported"], [
    AC_MSG_WARN([Disabling OpenMP support since no suitable C compiler found.])
])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset])

AC_CONFIG_FILES([Makefile
                 examples/Makefile])
AC_OUTPUT
