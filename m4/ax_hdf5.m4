#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_HDF5],[
	AC_ARG_WITH([hdf5],
		    AS_HELP_STRING([--with-hdf5],
		    [Compile with support for HDF5]),
		    [with_hdf5=${withval}], [with_hdf5=no])
                    
        if test "x${with_hdf5}" != xno; then
           PKG_CHECK_MODULES([HDF5_Fortran],[hdf5_fortran >= 1.14.0],
                             have_hdf5=yes, have_hdf5=no)
	   if test "x${have_hdf5}" = xyes; then 	
              LIBS="$HDF5_Fortran_LIBS $LIBS"
              FCFLAGS="$HDF5_Fortran_CFLAGS $FCFLAGS"
              AC_DEFINE(HAVE_HDF5,[1],
                     [Define if you have the HDF5 library.])
           fi
	fi
])
