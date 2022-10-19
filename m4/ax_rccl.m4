#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_RCCL],[
        AC_ARG_WITH([rccl],
                    AS_HELP_STRING([--with-rccl=DIR],
                    [Compile with support for RCCL]),
                    [
                    if test -d "$withval"; then
                       ac_rccl_path="$withval";
                       RCCL_LDFLAGS="-L$ac_rccl_path/lib"
                       RCCL_CFLAGS="-I$ac_rccl_path/include"
                    fi
                    ], [with_rccl=no])
        rccl_bcknd="0"
        if test "x${with_rccl}" != xno; then
           if test "x${have_hip}" != xyes; then
              AC_MSG_ERROR([RCCL requires HIP])
           fi
           
           if test -d "$ac_rccl_path"; then
	      CPPFLAGS_SAVED="$CPPFLAGS"
              LDFLAGS_SAVED="$LDFLAGS"
	      CPPFLAGS="$RCCL_CPPFLAGS $CPPFLAGS"
	      LDFLAGS="$RCCL_LDFLAGS $LDFLAGS"
	      export CPPFLAGS
	      export LDFLAGS
	   fi

              _CC=$CC
	      _LIBS=$LIBS
	      AC_LANG_PUSH([C])
	      CC=$HIPCC
	      LIBS=""

	      AC_CHECK_LIB(rccl, ncclCommDestroy,
	                   [have_rccl=yes;RCCL_LIBS="-lrccl"],
                           [have_rccl=no],[])
              AC_SUBST(have_rccl)
              if test x"${have_rccl}" = xyes; then
                 rccl_bcknd="1"
                 AC_DEFINE(HAVE_RCCL,1,[Define if you have RCCL.])
                 LIBS="$RCCL_LIBS $_LIBS"
                 HIP_HIPCC_FLAGS="$HIP_HIPCC_FLAGS $RCCL_CFLAGS"
              else
                 AC_MSG_ERROR([RCCL not found])
              fi
              CC=$_CC
              AC_LANG_POP([C])
        fi
        AC_SUBST(rccl_bcknd)      
])
