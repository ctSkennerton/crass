#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(crass, 0.0.9, c.skennerton@gmail.com)
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/crass.cpp])
AC_CONFIG_HEADER([config.h])

# custom defines for config.h
AC_DEFINE([PACKAGE_FULL_NAME],["CRisprASSembler"],[Define the full name of the package])
AC_DEFINE([PACKAGE_MAJOR_VERSION],[0],[Define the major version of the package])
AC_DEFINE([PACKAGE_MINOR_VERSION],[0],[Define the minor version of the package])
AC_DEFINE([PACKAGE_REVISION],[9],[Define the revision of the package])

# extra configure options
#
#

#DEBUG
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
    [Turn on debugging preprocessor flags])],
    [case "${enableval}" in
        yes) debug=true ;;
        no)  debug=false ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
    esac],[debug=false])

AM_CONDITIONAL([DEBUG],[test x$debug = xtrue])

if test x$debug = xtrue; then
    AC_DEFINE([DEBUG],[1],[define to 1 for debug statments])
fi

#SUPER_LOGGING
AC_ARG_ENABLE([verbose_logging],
    [AS_HELP_STRING([--enable-verbose_logging],
    [Change the verbosity of log file statments to include the function and line numbers])],
    [case "${enableval}" in
        yes) verbose_logging=true ;;
        no)  verbose_logging=false ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-verbose_logging]) ;;
    esac],[verbose_logging=false])

if test x$verbose_logging = xtrue; then
    AC_DEFINE([SUPER_LOGGING],[1],[define to 1 for verbose logging statments])
fi

# enable rendering of images with graphviz 
AC_ARG_ENABLE([rendering],
    [AS_HELP_STRING([--enable-rendering],
    [Enables images being drawn in the output files])],
    [case "${enableval}" in
        yes) rendering=true ;;
        no)  rendering=false ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-rendering]) ;;
    esac],[rendering=false])

#check that graphviz is installed before defining RENDERING
if test x$rendering = xtrue; then
    AC_CHECK_PROGS([GV_PROGS], [neato dot circo fdp sfdp twopi], [no])
    if test x$GV_PROGS != xno; then
        AC_DEFINE([RENDERING],[1],[Defines to 1 if image rendering is to be performed])
    fi
fi



# Checks for different libs, headers, functions & programs
#
#

# Checks for compilers.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

AX_XERCESC

AC_CHECK_LIB([z],[gzread],,[AC_MSG_ERROR([zlib not found])])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h stdlib.h string.h unistd.h getopt.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow strchr strspn])

# Checks for Graphviz program neato
AC_CHECK_PROG([NEATO],[neato],[yes],[no])
if test x$NEATO = xyes; then
    AC_DEFINE([HAVE_NEATO],[1],[Defines to 1 if the program 'neato' is found in PATH])
fi 

# Checks for Graphviz program dot
AC_CHECK_PROG([DOT],[dot],[yes],[no])
if test x$DOT = xyes; then
    AC_DEFINE([HAVE_DOT],[1],[Defines to 1 if the program 'dot' is found in PATH])
fi 

# Checks for Graphviz program fdp
AC_CHECK_PROG([FDP],[fdp],[yes],[no])
if test x$FDP = xyes; then
    AC_DEFINE([HAVE_FDP],[1],[Defines to 1 if the program 'fdp' is found in PATH])
fi 

# Checks for Graphviz program sfdp
AC_CHECK_PROG([SFDP],[sfdp],[yes],[no])
if test x$SFDP = xyes; then
    AC_DEFINE([HAVE_SFDP],[1],[Defines to 1 if the program 'sfdp' is found in PATH])
fi 

# Checks for Graphviz program circo
AC_CHECK_PROG([CIRCO],[circo],[yes],[no])
if test x$CIRCO = xyes; then
    AC_DEFINE([HAVE_CIRCO],[1],[Defines to 1 if the program 'circo' is found in PATH])
fi 

# Checks for Graphviz program twopi
AC_CHECK_PROG([TWOPI],[twopi],[yes],[no])
if test x$FDP = xyes; then
    AC_DEFINE([HAVE_TWOPI],[1],[Defines to 1 if the program 'twopi' is found in PATH])
fi 


# Checks for Cap3 assembler
AC_CHECK_PROG([CAP],[cap3],[yes],[no])
if test x$CAP = xyes; then
    AC_DEFINE([HAVE_CAP3],[1],[Defines to 1 if the program 'cap3' is found in PATH])
fi 

# Checks for both executables for velvet
AC_CHECK_PROG([VELVET_G],[velvetg],[yes],[no])
AC_CHECK_PROG([VELVET_H],[velveth],[yes],[no])

if test x$VELVET_G = xyes; then
    if test x$VELVET_H = xyes; then
        AC_DEFINE([HAVE_VELVET],[1],[Defines to 1 if the program 'velvetg' and 'velveth' is found in PATH])
    fi
fi 

AC_OUTPUT(Makefile src/Makefile man/Makefile)
