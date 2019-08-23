#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([argot-io], [1.0.0])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([README.md])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
AX_MPI

# Checks for libraries.
AC_CHECK_LIB([m], [exp])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT64_T
AC_TYPE_MODE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([getcwd gettimeofday mkdir pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
