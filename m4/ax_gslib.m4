#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_GSLIB],[
	AC_ARG_WITH([gslib],
                    AS_HELP_STRING([--with-gslib=DIR],
                    [Compile with gslib]),
                    [	   
                    if test -d "$withval"; then
                       ac_gslib_path="$withval";
               	       AS_IF([test -d "$ac_cuda_path/lib64"],
          	             [suffix="64"],[suffix=""])
		       GSLIB_LDFLAGS="-L$ac_gslib_path/lib$suffix -L$ac_gslib_path"
                    fi
                    ], [with_gslib=no])
        if test "x${with_gslib}" != xno; then
            	if test -d "$ac_gslib_path"; then
                   LDFLAGS_SAVED="$LDFLAGS"
                   LDFLAGS="$GSLIB_LDFLAGS $LDFLAGS"
                   export LDFLAGS
                fi

                _LIBS=$LIBS
                LIBS=""

                AC_CHECK_LIB(gs, fgslib_gs_setup,
		             [have_gslib=yes;GSLIB_LIBS="-lgs"],
		             [have_gslib=no])
                AC_SUBST(have_gslib)
                if test x"${have_gslib}" = xyes; then
                   AC_DEFINE(HAVE_GSLIB,1,[Define if you have the GS library.])
                   LIBS="$GSLIB_LIBS $_LIBS"
                else
                   AC_MSG_ERROR([gslib not found])
		fi
	fi
])



