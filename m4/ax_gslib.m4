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
	[Directory for gslib]),
	[	   
	if test -d "$withval"; then
		ac_gslib_path="$withval";
		GSLIB_LDFLAGS="-L$ac_gslib_path/"  
	fi
	],)

	if test -d "$ac_gslib_path"; then
	   LDFLAGS_SAVED="$LDFLAGS"
	   LDFLAGS="$GSLIB_LDFLAGS $LDFLAGS"
	   export LDFLAGS
	fi

	AC_LANG(Fortran)

	_LIBS=$LIBS

	AC_CHECK_LIB(gs, fgslib_gs_setup,
			       [have_gslib=yes;GSLIB_LIBS="-lgs"],
			       [have_gslib=no])
	AC_SUBST(GSLIB_LIBS)
	if test x"${have_gslib}" = xyes; then
	   AC_DEFINE(HAVE_GSLIB,1,[Define if you have the GS library.])
	   LIBS="$GSLIB_LIBS $_LIBS"
	else
		if test -d "$ac_gslib_path"; then	
		   LDFLAGS="$LDFLAGS_SAVED"
		fi
	fi

])



