#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(crass, 1.0.0-alpha, c.skennerton@gmail.com)

AC_CONFIG_AUX_DIR(build)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/crass/crass.cpp])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([foreign])
LT_INIT
LT_LANG([C++])


# Checks for different libs, headers, functions & programs
#
#

# Checks for compilers.


AC_PROG_CXX([clang++ g++ c++])

AC_CANONICAL_HOST
case $host_os in
    darwin* )
    CXX=clang++
    CXXFLAGS+=" -stdlib=libc++"
    ;;
esac
#AX_LIBCRISPR
#if test $HAVE_LIBCRISPR = no; then
#AC_MSG_ERROR([Cannot find licrispr])
#else
#AM_CONDITIONAL([LIBCRISPR_FOUND], [test $HAVE_LIBCRISPR = yes])
#fi

PERFORM_CRASS_ASSEMBLY=no
# Checks for libraries.
#if test $assembly = true; then
#    AC_CHECK_PROGS([ASSEMBLERS], [velvetg cap3], [no])
#    if test ASSEMBLERS != no; then
#        PERFORM_CRASS_ASSEMBLY=yes
#        AC_DEFINE([PERFORM_CRASS_ASSEMBLY],[1],[Defines to 1 if your system configuration enables you to use the Crass assembler])
#    else
#        echo "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"
#        echo "No genome assemblers were found"
#        echo "in your PATH!"
#        echo "Your installation will not conatin"
#        echo "the Crass assembler!"
#        echo "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"
#    fi
#fi
#AM_CONDITIONAL([ASSEMBLY_WRAPPER], [test x$PERFORM_CRASS_ASSEMBLY = xyes])

zlib_flags=
AC_CHECK_LIB([z],[gzread],
             [zlib_flags="-lz"
              AC_DEFINE([HAVE_ZLIB],[1],[Defines to 1 if the z library (-lz) is found])],
             [AC_MSG_ERROR([zlib not found])])
AC_SUBST(zlib_flags)
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([getopt.h])




# Check for presence of pdfLaTeX
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
if test -z "$PDFLATEX"; then
AC_MSG_WARN([Unable to create PDF version of the user manual.])
fi

AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")


# only perform the check if rendering
if test x$rendering = xtrue; then
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
fi

AC_OUTPUT(Makefile src/Makefile src/crass/Makefile man/Makefile doc/Makefile)
