#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_CUDA],[
	AC_ARG_WITH([cuda],
		    AS_HELP_STRING([--with-cuda=DIR],
		    [Compile with CUDA backend]),
		    [
		    if test -d "$withval"; then
		       ac_cuda_path="$withval";
       		       AS_IF([test -d "$ac_cuda_path/lib64"],
          	             [suffix="64"],[suffix=""])
		       CUDA_LDFLAGS="-L$ac_cuda_path/lib$suffix"
		    fi
		    ], [with_cuda=no])
	cuda_bcknd="0"
	if test "x${with_cuda}" != xno; then
	        if test -d "$ac_cuda_path"; then
	   	   CPPFLAGS_SAVED="$CPPFLAGS"
		   LDFLAGS_SAVED="$LDFLAGS"
		   CPPFLAGS="$CUDA_CPPFLAGS $CPPFLAGS"
		   LDFLAGS="$CUDA_LDFLAGS"
		   export CPPFLAGS
		   export LDFLAGS
		   AC_PATH_PROG(NVCC, nvcc, "no")
		fi
		
                AS_IF([test "$CUDA_CFLAGS"],[],[CUDA_CFLAGS="-O3"])
		
		_CC=$CC
		_LIBS=$LIBS
		AC_LANG_PUSH([C])
		CC=$NVCC
		LIBS=""

		AC_CHECK_LIB(cudart, cudaFree,
		             [have_cuda=yes;CUDA_LIBS="-lcudart"],[have_cuda=no])
		AC_SUBST(have_cuda)		
		if test x"${have_cuda}" = xyes; then		   
                   cuda_bcknd="1"
		   AC_DEFINE(HAVE_CUDA,1,[Define if you have CUDA.])
		   LIBS="$CUDA_LIBS $_LIBS"
		   LDFLAGS="$CUDA_LDFLAGS $LDFLAGS_SAVED"
		else
		   AC_MSG_ERROR([CUDA not found])
		fi
	        CC=$_CC
		AC_LANG_POP([C])
	fi
	AC_SUBST(cuda_bcknd)	
])
