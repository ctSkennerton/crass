#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(crass, 0.3.12, c.skennerton@gmail.com)

AC_CONFIG_AUX_DIR(build)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/crass/crass.cpp])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([foreign])
LT_INIT
LT_LANG([C++])
# custom defines for config.h
AC_DEFINE([PACKAGE_FULL_NAME],["CRisprASSembler"],[Define the full name of the package])
AC_DEFINE([PACKAGE_MAJOR_VERSION],[0],[Define the major version of the package])
AC_DEFINE([PACKAGE_MINOR_VERSION],[3],[Define the minor version of the package])
AC_DEFINE([PACKAGE_REVISION],[12],[Define the revision of the package])


# extra configure options
#
#

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

# check for assembly
AC_ARG_ENABLE([assembly],
    [AS_HELP_STRING([--disable-assembly],
        [Enables images being drawn in the output files])],
        [case "${enableval}" in
            yes) assembly=true ;;
            no)  assembly=false ;;
            *) AC_MSG_ERROR([bad value ${enableval} for --disable-assembly]) ;;
        esac],[assembly=true])

#check that graphviz is installed before defining RENDERING

if test x$rendering = xtrue; then
    AC_CHECK_PROGS([GV_PROGS], [neato fdp sfdp dot circo twopi], [no])

    if test x$GV_PROGS != xno; then
        
        AC_DEFINE([RENDERING],[1],[Defines to 1 if image rendering is to be performed])
        AC_DEFINE_UNQUOTED([DEFAULT_RENDERING_ALGORITHM],["${GV_PROGS}"],[Sets the Default Graphviz algorithm])
    else
        echo "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"
        echo "None of the Graphviz executables could"
        echo "be found in your PATH!"
        echo "Your installation will not conatin"
        echo "rendering options!"
        echo "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"
    fi
fi


# Checks for different libs, headers, functions & programs
#
#

# Checks for compilers.


AC_PROG_CXX
AC_PROG_CC

AX_LIB_XERCES

if test $HAVE_XERCES = no; then
AC_MSG_ERROR([Cannot find Xerces-c])
else
AM_CONDITIONAL([XERCES_FOUND], [test $HAVE_XERCES = yes])
fi

PERFORM_CRASS_ASSEMBLY=no 
# Checks for libraries.
if test $assembly = true; then
    AC_CHECK_PROGS([ASSEMBLERS], [velvetg cap3], [no])
    if test ASSEMBLERS != no; then
        PERFORM_CRASS_ASSEMBLY=yes
        AC_DEFINE([PERFORM_CRASS_ASSEMBLY],[1],[Defines to 1 if your system configuration enables you to use the Crass assembler])
    else
        echo "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"
        echo "No genome assemblers were found"
        echo "in your PATH!"
        echo "Your installation will not conatin"
        echo "the Crass assembler!"
        echo "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"
    fi
fi
AM_CONDITIONAL([ASSEMBLY_WRAPPER], [test x$PERFORM_CRASS_ASSEMBLY = xyes])

zlib_flags=
AC_CHECK_LIB([z],[gzread],
             [zlib_flags="-lz"
              AC_DEFINE([HAVE_ZLIB],[1],[Defines to 1 if the z library (-lz) is found])],
             [AC_MSG_ERROR([zlib not found])])
AC_SUBST(zlib_flags)
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
#check that graphviz is installed before defining RENDERING
GV_LIBS=
GV_BOOL=false
if test x$rendering = xtrue; then
    AC_CHECK_PROGS([GV_PROGS], [neato fdp sfdp dot circo twopi], [no])
    #check that graphviz is installed before defining RENDERING
    AC_CHECK_LIB([gvc],[gvContext],
                 [GV_LIBS="-lgvc $GV_LIBS" 
                  GV_BOOL=true
                 AC_DEFINE([HAVE_LIBGVC],[1],[Defines to one if the gvc library (-lgvc) is found])],
                 [GV_BOOL=false])
    AC_CHECK_LIB([cdt],[dtopen],
                 [GV_LIBS="-lcdt $GV_LIBS"
                  GV_BOOL=true
                  AC_DEFINE([HAVE_LIBCDT],[1],[Defines to one if the cdt library (-lcdt) is found])],
                 [GV_BOOL=false])
    AC_CHECK_LIB([graph],[agnode],
                 [GV_LIBS="-lgraph $GV_LIBS"
                  GV_BOOL=true
                  AC_DEFINE([HAVE_LIBGRAPH],[1],[Defines to one if the graph library (-lgraph) is found])],
                 [GV_BOOL=false])
    if test x$GV_PROGS != xno; then
        
        AC_DEFINE([RENDERING],[1],[Defines to 1 if image rendering is to be performed])
        AC_DEFINE_UNQUOTED([DEFAULT_RENDERING_ALGORITHM],["${GV_PROGS}"],[Sets the Default Graphviz algorithm])
    else
        echo "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"
        echo "None of the Graphviz executables could"
        echo "be found in your PATH!"
        echo "Your installation will not conatin"
        echo "rendering options!"
        echo "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"
    fi
fi
AM_CONDITIONAL([FOUND_GRAPHVIZ_LIBRARIES],[test $GV_BOOL = true])
AC_SUBST(GV_LIBS)


if test $PERFORM_CRASS_ASSEMBLY = yes; then
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
fi

AC_OUTPUT(Makefile src/Makefile src/crass/Makefile man/Makefile doc/Makefile)
