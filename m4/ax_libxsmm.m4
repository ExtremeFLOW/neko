#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_LIBXSMM],[
	AC_ARG_WITH([libxsmm],
		    AS_HELP_STRING([--with-libxsmm],
		    [Compile with support for libxsmm]),
		    [with_libxsmm=${withval}], [with_libxsmm=no])

	xsmm_bcknd="0"
        if test "x${with_libxsmm}" != xno; then
     	   PKG_CHECK_MODULES([libxsmmf], [libxsmmf >= 0.16.1], 
	           	     have_libxsmm=yes xsmm_bcknd="1",
			     xsmm_bkcnd="0" have_libxsmm=no)
	   if test "x${have_libxsmm}" = xyes; then 	
	      FCFLAGS="$FCFLAGS $libxsmmf_CFLAGS"
	      LIBS="$LIBS $libxsmmf_LIBS"
	      AC_DEFINE(HAVE_LIBXSMM,[1],
			[Define if you have the LIBXSMM library.])
	   fi
	fi
        AC_SUBST(xsmm_bcknd)

])
