#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_NCCL],[
        AC_ARG_WITH([nccl],
                    AC_HELP_STRING([--with-nccl=DIR],
                    [Compile with support for NCCL]),
                    [
                    if test -d "$withval"; then
                       ac_nccl_path="$withval";
                       AS_IF([test -d "$ac_nccl_path/lib64"],
                             [suffix="64"],[suffix=""])
                       NCCL_LDFLAGS="-L$ac_nccl_path/lib$suffix"
                       NCCL_CFLAGS="-I$ac_nccl_path/include"
                    fi
                    ], [with_nccl=no])
        nccl_bcknd="0"
        if test "x${with_nccl}" != xno; then
           if test "x${have_cuda}" != xyes; then
              AC_MSG_ERROR([NCCL requires CUDA])
           fi
           
           if test -d "$ac_nccl_path"; then
	      CPPFLAGS_SAVED="$CPPFLAGS"
              LDFLAGS_SAVED="$LDFLAGS"
	      CPPFLAGS="$NCCL_CPPFLAGS $CPPFLAGS"
	      LDFLAGS="$NCCL_LDFLAGS $LDFLAGS"
	      export CPPFLAGS
	      export LDFLAGS
	   fi

              _CC=$CC
	      _LIBS=$LIBS
	      AC_LANG_PUSH([C])
	      CC=$NVCC
	      LIBS=""

	      AC_CHECK_LIB(nccl, ncclCommDestroy,
	                   [have_nccl=yes;NCCL_LIBS="-lnccl"],
                           [have_nccl=no],[])
              AC_SUBST(have_nccl)
              if test x"${have_nccl}" = xyes; then
                 nccl_bcknd="1"
                 AC_DEFINE(HAVE_NCCL,1,[Define if you have NCCL.])
                 LIBS="$NCCL_LIBS $_LIBS"
                 CUDA_CFLAGS="$CUDA_CFLAGS $NCCL_CFLAGS"
              else
                 AC_MSG_ERROR([NCCL not found])
              fi
              CC=$_CC
              AC_LANG_POP([C])
        fi
        AC_SUBST(nccl_bcknd)      
])
