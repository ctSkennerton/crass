dnl SYNOPSIS
dnl
dnl   AX_LIBCRISPR
dnl
dnl DESCRIPTION
dnl
dnl   This macro provides tests of availability for libcrispr. This macros 
dnl   checks for libcrispr Parser headers and libraries and defines compilation flags
dnl
dnl   Macro supports following options and their values:
dnl
dnl   1) Single-option usage:
dnl
dnl     --with-libcrispr - yes, no or path to libcripsr installation prefix
dnl
dnl   2) Three-options usage (all options are required):
dnl
dnl     --with-libcripr=yes
dnl     --with-libcrispr-inc - path to base directory with libcrispr headers
dnl     --with-libcrispr-lib - linker flags for libcrispr
dnl
dnl   This macro calls:
dnl
dnl     AC_SUBST(LIBCRISPR_CPPFLAGS)
dnl     AC_SUBST(LIBCRISPR_LDFLAGS)
dnl     AC_SUBST([LIBCRISPR_LIBS])
dnl
dnl   And sets:
dnl
dnl     HAVE_LIBCRISPR
dnl
dnl LICENSE
dnl
dnl   Copyright (c) 2012 Connor Skennerton
dnl
dnl   Copying and distribution of this file, with or without modification, are
dnl   permitted in any medium without royalty provided the copyright notice
dnl   and this notice are preserved. This file is offered as-is, without any
dnl   warranty.


AC_DEFUN([AX_LIBCRISPR],
[
    AC_REQUIRE([AX_LIB_XERCES])

    AC_ARG_WITH([libcrispr],
        AS_HELP_STRING([--with-libcrispr=@<:@ARG@:>@],
            [use libcrispr from given prefix (ARG=path); check standard prefixes (ARG=yes); disable (ARG=no)]
        ),
        [
        if test "$withval" = "yes"; then
            if test -d /usr/local/include/libcrispr ; then
                libcrispr_prefix=/usr/local
            elif test -d /usr/include/libcrispr ; then
                libcrispr_prefix=/usr
            else
                libcrispr_prefix=""
            fi
            libcrispr_requested="yes"
        elif test -d "$withval"; then
            libcrispr_prefix="$withval"
            libcrispr_requested="yes"
        else
            libcrispr_prefix=""
            libcrispr_requested="no"
        fi
        ],
        [
        dnl Default behavior is implicit yes
        if test -d /usr/local/include/libcrispr ; then
            libcrispr_prefix=/usr/local
        elif test -d /usr/include/libcrispr ; then
            libcrispr_prefix=/usr
        else
            libcrispr_prefix=""
        fi
        ]
    )

    AC_ARG_WITH([libcrispr-inc],
        AS_HELP_STRING([--with-libcrispr-inc=@<:@DIR@:>@],
            [path to libcrispr headers]
        ),
        [libcrispr_include_dir="$withval"],
        [libcrispr_include_dir=""]
    )
    AC_ARG_WITH([libcrispr-lib],
        AS_HELP_STRING([--with-libcrispr-lib=@<:@ARG@:>@],
            [link options for libcrispr libraries]
        ),
        [libcrispr_ldflags="$withval"],
        [libcrispr_ldflags=""]
    )

    LIBCRISPR_CPPFLAGS=""
    LIBCRISPR_LDFLAGS=""

    dnl
    dnl Collect include/lib paths and flags
    dnl
    run_libcrispr_test="no"

    if test -n "$libcrispr_prefix"; then
        libcrispr_include_dir="$libcrispr_prefix/include"
        libcrispr_include_dir2="$libcrispr_prefix/include/libcrispr"
        libcrispr_ldflags="-L$libcrispr_prefix/lib"
        run_libcrispr_test="yes"
    elif test "$libcrispr_requested" = "yes"; then
        if test -n "$libcrispr_include_dir" -a -n "$libcrispr_lib_flags"; then
            libcrispr_include_dir2="$libcrispr_include_dir/libcrispr"
            run_libcrispr_test="yes"
        fi
    else
        run_libcrispr_test="no"
    fi

    libcrispr_libs="-lcrispr"
    dnl
    dnl Check libcrispr files
    dnl

    if test "$run_libcrispr_test" = "yes"; then

        saved_CPPFLAGS="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS -I$libcrispr_include_dir -I$libcrispr_include_dir2 $XERCES_CPPFLAGS" 

        saved_LDFLAGS="$LDFLAGS"
        LDFLAGS="$LDFLAGS $libcrispr_ldflags $XERCES_LDFLAGS $PTHREAD_LDFLAGS"

        saved_LIBS="$LIBS"
        LIBS="$libcrispr_libs $XERCES_LIBS $PTHREAD_LIBS $LIBS"
        dnl echo "libs: ${LIBS}"
        dnl echo "cppflags: ${CPPFLAGS}"
        dnl echo "ldflags: ${LDFLAGS}"
        dnl
        dnl Check libcrispr headers
        dnl
        AC_MSG_CHECKING([for libcrispr headers in $libcrispr_include_dir and $libcrispr_include_dir2])

        AC_LANG_PUSH([C++])
        AC_COMPILE_IFELSE([
            AC_LANG_PROGRAM(
                [[
@%:@include <libcrispr/Exception.h>
@%:@include <libcrispr/base.h>
@%:@include <libcrispr/writer.h>
@%:@include <libcrispr/reader.h>
@%:@include <libcrispr/parser.h>
@%:@include <libcrispr/StlExt.h>
                ]],
                [[]]
            )],
            [
            LIBCRISPR_CPPFLAGS="-I$libcrispr_include_dir -I$libcrispr_include_dir2"
            libcrispr_header_found="yes"
            AC_MSG_RESULT([found])
            ],
            [
            libcrispr_header_found="no"
            AC_MSG_RESULT([not found])
            ]
        )
        AC_LANG_POP([C++])

        dnl
        dnl Check libcrispr libraries
        dnl
        if test "$libcrispr_header_found" = "yes"; then

            AC_MSG_CHECKING([for libcrispr libraries])

            AC_LANG_PUSH([C++])
            AC_LINK_IFELSE([
                AC_LANG_PROGRAM(
                    [[
@%:@include <libcrispr/base.h>
                    ]],
                    [[
crispr::xml::base crispr_test;
                    ]]
                )],
                [
                LIBCRISPR_LDFLAGS="$libcrispr_ldflags $XERCES_LDFLAGS $PTHREAD_LDFLAGS"
                LIBCRISPR_LIBS="$libcrispr_libs $XERCES_LIBS $PTHREAD_LIBS"
                libcrispr_lib_found="yes"
                AC_MSG_RESULT([found])
                ],
                [
                libcrispr_lib_found="no"
                AC_MSG_RESULT([not found])
                ]
            )
            AC_LANG_POP([C++])
        fi

        CPPFLAGS="$saved_CPPFLAGS"
        LDFLAGS="$saved_LDFLAGS"
        LIBS="$saved_LIBS"
    fi

    AC_MSG_CHECKING([for libcrispr])

    if test "$run_libcrispr_test" = "yes"; then
        if test "$libcrispr_header_found" = "yes" -a "$libcrispr_lib_found" = "yes"; then

            AC_SUBST([LIBCRISPR_CPPFLAGS])
            AC_SUBST([LIBCRISPR_LDFLAGS])
            AC_SUBST([LIBCRISPR_LIBS])

            HAVE_LIBCRISPR="yes"
        else
            HAVE_LIBCRISPR="no"
        fi

        AC_MSG_RESULT([$HAVE_LIBCRISPR])


    else
        HAVE_LIBCRISPR="no"
        AC_MSG_RESULT([$HAVE_LIBCRISPR])

        if test "$libcrispr_requested" = "yes"; then
            AC_MSG_WARN([libcrispr requested but headers or library not found. Specify valid prefix of libcrispr using --with-libcrispr=@<:@DIR@:>@ or provide include directory and linker flags using --with-libcrispr-inc and --with-libcrispr-lib])
        fi
    fi
])
