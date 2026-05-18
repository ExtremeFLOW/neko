#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_NVTX],[
        AC_ARG_WITH([nvtx],
                    AS_HELP_STRING([--with-nvtx=DIR],
                    [Compile with support for NVTX]),
                    [
                    if test -d "$withval"; then
                       ac_nvtx_path="$withval";
                       AS_IF([test -d "$ac_nvtx_path/lib64"],
                             [suffix="64"],[suffix=""])
                       NVTX_LDFLAGS="-L$ac_nvtx_path/lib$suffix"
                       NVTX_CFLAGS="-I$ac_nvtx_path/include"
                    fi
                    ], [with_nvtx=no])
        nvtx_bcknd="0"
        if test "x${with_nvtx}" != xno; then
           if test "x${have_cuda}" != xyes; then
              AC_MSG_ERROR([NVTX requires CUDA])
           fi
           
           if test -d "$ac_nvtx_path"; then
	      CPPFLAGS_SAVED="$CPPFLAGS"
              LDFLAGS_SAVED="$LDFLAGS"
	      CPPFLAGS="$NVTX_CPPFLAGS $CPPFLAGS"
	      LDFLAGS="$NVTX_LDFLAGS $LDFLAGS"
	      export CPPFLAGS
	      export LDFLAGS
	   fi

              _CC=$CC
              _CFLAGS=$CFLAGS
	      _LIBS=$LIBS
	      AC_LANG_PUSH([C])
	      CC=$NVCC
              CFLAGS=""
	      LIBS=""

	      AC_CHECK_LIB(nvToolsExt, nvtxRangePop,
	                   [have_nvtx=yes;NVTX_LIBS="-lnvToolsExt"],
                           [have_nvtx=no],[])

              # CUDA 12.9+ drops libnvToolsExt; fall back to header-only NVTX3
              if test x"${have_nvtx}" != xyes; then
                 AC_MSG_CHECKING([for NVTX3 header-only library])
                 AC_COMPILE_IFELSE(
                    [AC_LANG_PROGRAM([[#include <nvtx3/nvToolsExt.h>]],
                                     [[nvtxRangePop();]])],
                    [have_nvtx=yes; have_nvtx3=yes; NVTX_LIBS=""
                     AC_MSG_RESULT([yes])],
                    [AC_MSG_RESULT([no])])
              fi

              AC_SUBST(have_nvtx)
              if test x"${have_nvtx}" = xyes; then
                 nvtx_bcknd="1"
                 AC_DEFINE(HAVE_NVTX,1,[Define if you have NVTX.])
                 NVTX_CFLAGS="$NVTX_CFLAGS -DHAVE_NVTX"
                 if test x"${have_nvtx3}" = xyes; then
                    AC_DEFINE(HAVE_NVTX3,1,[Define if NVTX3 (header-only) is used.])
                    NVTX_CFLAGS="$NVTX_CFLAGS -DHAVE_NVTX3"
                 fi
                 LIBS="$NVTX_LIBS $_LIBS"
                 CUDA_CFLAGS="$CUDA_CFLAGS $NVTX_CFLAGS"
              else
                 AC_MSG_ERROR([NVTX not found])
              fi
              CC=$_CC
              CFLAGS=$_CFLAGS
              AC_LANG_POP([C])
        fi
        AC_SUBST(nvtx_bcknd)      
])
