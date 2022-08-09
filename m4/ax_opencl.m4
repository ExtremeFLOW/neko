#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_OPENCL],[
	AC_ARG_WITH([opencl],
		    AC_HELP_STRING([--with-opencl=DIR],
		    [Compile with OpenCL backend]),
		    [
   		    if test -d "$withval"; then
		       ac_opencl_path="$withval";
         	       AS_IF([test -d "$ac_opencl_path/lib64"],
          	             [suffix="64"],[suffix=""])
     		       OPENCL_LDFLAGS="-L$ac_opencl_path/lib$suffix"
       		       OPENCL_CPPFLAGS="-I$ac_opencl_path/include"
		       OPENCL_LIB="-lOpenCL"

		    else
		    # Assume we're on a Mac unless an OpenCL dir is given
		       OPENCL_LIB="-framework OpenCL"
		    fi
		    ], [with_opencl=no])
	opencl_bcknd="0"
	if test "x${with_opencl}" != xno; then
	        if test -d "$ac_opencl_path"; then
	   	   CPPFLAGS_SAVED="$CPPFLAGS"
		   LDFLAGS_SAVED="$LDFLAGS"
		   CPPFLAGS="$OPENCL_CPPFLAGS $CPPFLAGS"
		   LDFLAGS="$OPENCL_LDFLAGS $LDFLAGS"
		   export CPPFLAGS
		   export LDFLAGS
		fi

		AC_LANG_PUSH([C])
		AC_LANG_ASSERT([C])
		LIBS_SAVED="$LIBS"
		LIBS="$OPENCL_LIB $LIBS"
		AC_MSG_CHECKING([for OpenCL])
		AC_LINK_IFELSE([AC_LANG_SOURCE([
                #ifdef __APPLE__
                #include <OpenCL/opencl.h>
		#else
		#include <CL/cl.h>
		#endif
		int main(void) {
		    clGetPlatformIDs(0, NULL, NULL);
		}
		])],
		[have_opencl=yes],[have_opencl=no])
		
		AC_SUBST(have_opencl)		
		if test x"${have_opencl}" = xyes; then		   
                   opencl_bcknd="1"
		   AC_DEFINE(HAVE_OPENCL,1,[Define if you have OpenCL.])
		   AC_MSG_RESULT([yes])	
		else
		   AC_MSG_ERROR([OpenCL not found])
		fi
		AC_LANG_POP([C])
	fi
	AC_SUBST(opencl_bcknd)	
])
