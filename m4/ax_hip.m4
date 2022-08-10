#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_HIP],[
	AC_ARG_WITH([hip],
		    AC_HELP_STRING([--with-hip=DIR],
		    [Compile with HIP backend]),
		    [
		    if test -d "$withval"; then
		       ac_hip_path="$withval";
		       HIP_LDFLAGS="-L$ac_hip_path/lib"
		    fi
		    ], [with_hip=no])
	hip_bcknd="0"
	if test "x${with_hip}" != xno; then
	        if test -d "$ac_hip_path"; then
	   	   CPPFLAGS_SAVED="$CPPFLAGS"
		   LDFLAGS_SAVED="$LDFLAGS"
		   CPPFLAGS="$HIP_CPPFLAGS $CPPFLAGS"
		   LDFLAGS="$HIP_LDFLAGS $LDFLAGS"
		   export CPPFLAGS
		   export LDFLAGS
		   AS_IF([test "$HIPCC"],[],[AC_PATH_PROG(HIPCC, hipcc, "no")])
		fi

                AS_IF([test "$HIP_HIPCC_FLAGS"],[],[HIP_HIPCC_FLAGS="-O3"])

		_CC=$CC
		_LIBS=$LIBS
		AC_LANG_PUSH([C])
		CC=$HIPCC
		LIBS=""

		AC_CHECK_LIB(amdhip64, hipFree,
		             [have_hip=yes;HIP_LIBS="-lamdhip64"],[have_hip=no])
		AC_SUBST(have_hip)		
		if test x"${have_hip}" = xyes; then		   
                   hip_bcknd="1"
		   AC_DEFINE(HAVE_HIP,1,[Define if you have HIP.])
		   LIBS="$HIP_LIBS $_LIBS"
		else
		   AC_MSG_ERROR([HIP not found])
		fi
	        CC=$_CC
		AC_LANG_POP([C])
	fi
	AC_SUBST(hip_bcknd)	
])
