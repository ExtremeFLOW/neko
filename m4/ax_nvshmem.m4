#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_NVSHMEM],[
        AC_ARG_WITH([nvshmem],
                    AS_HELP_STRING([--with-nvshmem=DIR],
                    [Compile with support for NVSHMEM]),
                    [
                    if test -d "$withval"; then
                       ac_nvshmem_path="$withval";
                       AS_IF([test -d "$ac_nvshmem_path/lib64"],
                             [suffix="64"],[suffix=""])
                       NVSHMEM_LDFLAGS="-L$ac_nvshmem_path/lib$suffix"
                       NVSHMEM_CFLAGS="-I$ac_nvshmem_path/include"
                    fi
                    ], [with_nvshmem=no])
        nvshmem_bcknd="0"
        if test "x${with_nvshmem}" != xno; then
           if test "x${have_cuda}" != xyes; then
              AC_MSG_ERROR([NVSHMEM requires CUDA])
           fi
           
           if test -d "$ac_nvshmem_path"; then
	      CPPFLAGS_SAVED="$CPPFLAGS"
              LDFLAGS_SAVED="$LDFLAGS"
	      CPPFLAGS="$NVSHMEM_CPPFLAGS $CPPFLAGS"
	      LDFLAGS="$NVSHMEM_LDFLAGS $LDFLAGS"
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

	      AC_CHECK_LIB(nvshmem, nvshmem_quiet,
	                   [have_nvshmem=yes;NVSHMEM_LIBS="-lnvshmem"],
                           [have_nvshmem=no],[])
              AC_SUBST(have_nvshmem)
              if test x"${have_nvshmem}" = xyes; then
                 nvshmem_bcknd="1"
                 AC_DEFINE(HAVE_NVSHMEM,1,[Define if you have NVSHMEM.])
                 LIBS="$NVSHMEM_LIBS $_LIBS"
                 CUDA_CFLAGS="$CUDA_CFLAGS $NVSHMEM_CFLAGS -DHAVE_NVSHMEM=1"
                 CFLAGS="$_CFLAGS $NVSHMEM_CFLAGS"
              else
                 AC_MSG_ERROR([NVSHMEM not found])
                 CFLAGS=$_CFLAGS
              fi
              CC=$_CC

              AC_LANG_POP([C])
        fi
        AC_SUBST(nvshmem_bcknd)      
])
