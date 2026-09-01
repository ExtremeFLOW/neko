#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_UTOFU],[
        AC_ARG_WITH([utofu],
                    AS_HELP_STRING([--with-utofu=DIR],
                    [Compile with support for the native Tofu interconnect (uTofu)]),
                    [
                    if test -d "$withval"; then
                       ac_utofu_path="$withval";
                       UTOFU_LDFLAGS="-L$ac_utofu_path/lib"
                       UTOFU_CFLAGS="-I$ac_utofu_path/include"
                    fi
                    ], [with_utofu=no])
        have_utofu=no
        if test "x${with_utofu}" != xno; then
           CPPFLAGS_SAVED="$CPPFLAGS"
           LDFLAGS_SAVED="$LDFLAGS"
           CPPFLAGS="$UTOFU_CFLAGS $CPPFLAGS"
           LDFLAGS="$UTOFU_LDFLAGS $LDFLAGS"
           export CPPFLAGS
           export LDFLAGS

           AC_LANG_PUSH([C])
           AC_CHECK_HEADER([utofu.h],
                           [have_utofu_hdr=yes], [have_utofu_hdr=no])
           AC_CHECK_LIB([tofucom], [utofu_get_onesided_tnis],
                        [have_utofu=yes; UTOFU_LIBS="-ltofucom"],
                        [have_utofu=no], [])
           AC_LANG_POP([C])

           if test x"${have_utofu}" = xyes -a x"${have_utofu_hdr}" = xyes; then
              AC_DEFINE(HAVE_UTOFU, 1, [Define if you have uTofu.])
              LIBS="$UTOFU_LIBS $LIBS"
              CFLAGS="$CFLAGS $UTOFU_CFLAGS"
           else
              AC_MSG_ERROR([uTofu requested but utofu.h / libtofucom not found])
           fi
        fi
        AC_SUBST(have_utofu)
])
