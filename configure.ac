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
AC_OPENMP
AC_PROG_LIBTOOL
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

# Gfarm support
AC_MSG_CHECKING([for Gfarm file system])
AC_ARG_WITH([gfarm],
  [AS_HELP_STRING([--with-gfarm=GFARM_ROOT],
    [support IO with Gfarm backend @<:@default=no@:>@])],
  [], [with_gfarm=no])
AC_MSG_RESULT([$with_gfarm])
AM_CONDITIONAL([GFARM], [test x$with_gfarm != xno])
if test x$with_gfarm != xno; then
  case x$with_gfarm in
    xyes) ;;
    *)
    CPPFLAGS="$CPPFLAGS -I$with_gfarm/include"
    LDFLAGS="$LDFLAGS -L$with_gfarm/lib" ;;
  esac
  AC_CHECK_LIB([gfarm], [gfarm_initialize],, [AC_MSG_ERROR([libgfarm not found])])
fi

# CHFS support
PKG_CHECK_MODULES([CHFS], [chfs],
  [AC_DEFINE([HAVE_CHFS], [], [Build CHFS backend])
   CHFS_RPATH=$(pkg-config --libs-only-L chfs | sed 's/-L/-Wl,-rpath=/g')
   AC_SUBST(CHFS_RPATH)],
  [with_chfs=no])
AM_CONDITIONAL([CHFS], [test x$with_chfs != xno])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
