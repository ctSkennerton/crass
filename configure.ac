#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(crass, 0.0.1, c.skennerton@gmail.com)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/bm.cpp])
AC_CONFIG_HEADER([config.h])

# custom defines for config.h
AC_DEFINE([PACKAGE_FULL_NAME],["CRisprASSembler"],[Define the full name of the package])
AC_DEFINE([PACKAGE_MAJOR_VERSION],[0],[Define the major version of the package])
AC_DEFINE([PACKAGE_MINOR_VERSION],[0],[Define the minor version of the package])
AC_DEFINE([PACKAGE_REVISION],[1],[Define the revision of the package])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([z],[gzread],,[AC_MSG_ERROR([zlib not found in $LD_PATH])])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow strchr strspn])

AC_OUTPUT(Makefile src/Makefile man/Makefile)
