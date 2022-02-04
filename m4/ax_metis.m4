#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_METIS],[
	AC_ARG_WITH([metis],
	AS_HELP_STRING([--with-metis=DIR],
	[Directory for metis]),
	[	   
	if test -d "$withval"; then
		ac_metis_path="$withval";
		METIS_LDFLAGS="-L$ac_metis_path/lib"  
		METIS_CPPFLAGS="-I$ac_metis_path/include"
	fi
	],[with_metis=no])

	AC_ARG_WITH([metis-libdir],
	AS_HELP_STRING([--with-metis-libdir=LIBDIR],
	[Directory for metis library]),
	[
	if test -d "$withval"; then
	   ac_metis_libdir="$withval"
	fi
	],)

	if test -d "$ac_metis_libdir"; then	   
	    METIS_LDFLAGS="-L$ac_metis_libdir"  	   
        fi

	if test -d "$ac_metis_path"; then
	   CPPFLAGS_SAVED="$CPPFLAGS"
	   LDFLAGS_SAVED="$LDFLAGS"
	   CPPFLAGS="$METIS_CPPFLAGS $CPPFLAGS"
	   LDFLAGS="$METIS_LDFLAGS $LDFLAGS"
	   export CPPFLAGS
	   export LDFLAGS
	fi
	if test "x${with_metis}" != xno; then			
	   AC_LANG(C)
	   AC_CHECK_HEADER([metis.h],[have_metis_h=yes],[have_metis_h=no])
	   if test x"${have_metis_h}" = xno; then
		if test -d "$ac_metis_path"; then	
		   CPPFLAGS="$CPPFLAGS_SAVED"
		fi
	   fi
	   AC_LANG(Fortran)

	   AC_CHECK_LIB(metis, METIS_PartGraphKway,
			[have_metis=yes;METIS_LIBS="-lmetis"],
			[have_metis=no],[-lmetis])
	   AC_SUBST(METIS_LIBS)
	   if test x"${have_metis}" = xyes; then
	      AC_DEFINE(HAVE_METIS,1,[Define if you have the Metis library.])
	   else
	      if test -d "$ac_metis_path"; then	
	      	 LDFLAGS="$LDFLAGS_SAVED"
	      fi
	   fi
	fi

])



