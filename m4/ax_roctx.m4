#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_ROCTX],[
        AC_ARG_WITH([roctx],
                    AC_HELP_STRING([--with-roctx=DIR],
                    [Compile with support for ROCTX]),
                    [
                    if test -d "$withval"; then
                       ac_roctx_path="$withval";
                       ROCTX_LDFLAGS="-L$ac_roctx_path/lib"
                       ROCTX_CFLAGS="-I$ac_roctx_path/include"
                    fi
                    ], [with_roctx=no])
        roctx_bcknd="0"
        if test "x${with_roctx}" != xno; then
           if test "x${have_hip}" != xyes; then
              AC_MSG_ERROR([ROCTX requires HIP])
           fi
           
           if test -d "$ac_roctx_path"; then
	      CPPFLAGS_SAVED="$CPPFLAGS"
              LDFLAGS_SAVED="$LDFLAGS"
	      CPPFLAGS="$ROCTX_CPPFLAGS $CPPFLAGS"
	      LDFLAGS="$ROCTX_LDFLAGS $LDFLAGS"
	      export CPPFLAGS
	      export LDFLAGS
	   fi

              _CC=$CC
	      _LIBS=$LIBS
	      AC_LANG_PUSH([C])
	      CC=$HIPCC
	      LIBS=""

	      AC_CHECK_LIB(roctx64, roctxRangePop,
	                   [have_roctx=yes;ROCTX_LIBS="-lroctx64"],
                           [have_roctx=no],[])
              AC_SUBST(have_roctx)
              if test x"${have_roctx}" = xyes; then
                 roctx_bcknd="1"
                 AC_DEFINE(HAVE_ROCTX,1,[Define if you have ROCTX.])
                 LIBS="$ROCTX_LIBS $_LIBS"
                 CUDA_CFLAGS="$CUDA_CFLAGS $ROCTX_CFLAGS"
              else
                 AC_MSG_ERROR([ROCTX not found])
              fi
              CC=$_CC
              AC_LANG_POP([C])
        fi
        AC_SUBST(roctx_bcknd)      
])
