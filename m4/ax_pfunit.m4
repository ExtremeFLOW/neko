#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_PFUNIT], [
		      AC_ARG_WITH([pfunit],
		      AS_HELP_STRING([--with-pfunit=DIR],
		      [Directory for pFUnit]),
		      [
			if test -d "$withval"; then
			   PFUNIT_DIR="$withval";
			else
			   PFUNIT_DIR="/usr";
			fi
		      ],)
		      
 AC_REQUIRE([AC_PROG_FC])
 AC_LANG_PUSH([Fortran])


 AC_CHECK_FILE([$PFUNIT_DIR/include/PFUNIT.mk],
	       [have_pfunit=yes], [have_pfunit=no])

 AC_SUBST(PFUNIT_DIR)

 AM_CONDITIONAL([ENABLE_PFUNIT],[test "x${have_pfunit}" = xyes])
])
