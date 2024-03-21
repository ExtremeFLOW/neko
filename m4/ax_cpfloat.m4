#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_CPFLOAT],[
	AC_ARG_WITH([cpfloat],
                    AS_HELP_STRING([--with-cpfloat=DIR],
                    [Compile with cpfloat]),
                    [	   
                    if test -d "$withval"; then
                       ac_cpfloat_path="$withval";
               	       AS_IF([test -d "$ac_cpfloat_path/build/lib"],
          	             [suffix=""])
		       CPFLOAT_LDFLAGS="-L$ac_cpfloat_path/build/lib$suffix -L$ac_cpfloat_path/deps/pcg-c/src/"
                    fi
                    ], [with_cpfloat=no])
        if test "x${with_cpfloat}" != xno; then
            	if test -d "$ac_cpfloat_path"; then
		           LDFLAGS_SAVED="$LDFLAGS"
		           LDFLAGS="$LDFLAGS_SAVED $CPFLOAT_LDFLAGS"
		           export LDFLAGS
                fi

		        AC_LANG_PUSH([C])
		        CC=$CC

                AC_CHECK_LIB(cpfloat, init_optstruct,
		             [have_cpfloat=yes;CPFLOAT_LIBS="-lpcg_random -lcpfloat"],
		             [have_cpfloat=no])
                AC_SUBST(have_cpfloat)
                if test x"${have_cpfloat}" = xyes; then
                   AC_DEFINE(HAVE_CPFLOAT,1,[Define if you have the GS library.])
                   LIBS="$CPFLOAT_LIBS $LIBS"
                else
                   AC_MSG_ERROR([cpfloat not found])
		fi
	fi
])



