#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_MAMBA],[
	AC_ARG_WITH([mamba],
		    AC_HELP_STRING([--with-mamba],
		    [Compile with support for MAMBA]),
		    [with_mamba=${withval}], [with_mamba=no])


        if test "x${with_mamba}" != xno; then
     	   PKG_CHECK_MODULES([mamba], [mamba >= 0.1.6], 
	           	     have_mamba=yes, have_mamba=no)
           if test "x${have_mamba}" = xyes; then 	
	      FCFLAGS="$FCFLAGS $mamba_CFLAGS"
	      LIBS="$LIBS $mamba_LIBS -lfmmb"
	      AC_DEFINE(HAVE_MAMBA,[1],
			[Define if you have the MAMBA library.])
	   fi
	fi
])
