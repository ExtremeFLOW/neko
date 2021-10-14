#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_NOA],[
	AC_ARG_WITH([noa],
	AS_HELP_STRING([--with-noa=DIR],
	[Directory for Noa]),
	[	   
	if test -d "$withval"; then
		ac_noa_path="$withval";
		NOA_LDFLAGS="-L$ac_noa_path/"
		NOA_CPPFLAGS="-I$ac_noa_path/"  
	fi
	],)

	if test -d "$ac_noa_path"; then
	   CPPFLAGS_SAVED="$CPPFLAGS"
	   LDFLAGS_SAVED="$LDFLAGS"
   	   CPPFLAGS="$NOA_CPPFLAGS $CPPFLAGS"
	   LDFLAGS="$NOA_LDFLAGS $LDFLAGS"
	   export CPPFLAGS
	   export LDFLAGS	   
	fi
			

	AC_LANG_PUSH([C])

	PKG_CHECK_MODULES([PROTOBUF_C], [libprotobuf-c >= 1.0.0])
	CPPFLAGS="$PROTOBUF_C_CFLAGS $CPPFLAGS"
	LDFLAGS="$PROTOBUF_C_LDFLAGS $LDFLAGS"

	AC_CHECK_HEADER([noa.h],[have_noa_h=yes],[have_noa_h=no])
	if test x"${have_noa_h}" = xno; then
		if test -d "$ac_noa_path"; then	
		   CPPFLAGS="$CPPFLAGS_SAVED"
		fi
	fi

	AC_CHECK_LIB(noa, noa_init,
			       [have_noa=yes;NOA_LIBS="-lnoa"],
			       [have_noa=no])
	AC_SUBST(NOA_LIBS)
	if test x"${have_noa}" = xyes; then
	   AC_DEFINE(HAVE_NOA,1,[Define if you have the NoaSci library.])
	else
		if test -d "$ac_noa_path"; then	
		   LDFLAGS="$LDFLAGS_SAVED"
		fi
	fi
	AC_LANG_POP([C])
])



