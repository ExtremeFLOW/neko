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
           neko_hdf5_pc=yes
           if test -n "${HDF5_Fortran_LIBS}${HDF5_Fortran_CFLAGS}"; then
              neko_hdf5_pc=no
           fi
           PKG_CHECK_MODULES([HDF5_Fortran],[hdf5_fortran >= 1.14.0],
                             have_hdf5=yes, have_hdf5=no)
	   if test "x${have_hdf5}" = xyes; then
              LIBS="$HDF5_Fortran_LIBS $LIBS"
              NEKO_PKG_FCFLAGS="$HDF5_Fortran_CFLAGS $NEKO_PKG_FCFLAGS"
              if test "x${neko_hdf5_pc}" = xyes; then
                 NEKO_PKG_REQUIRES="$NEKO_PKG_REQUIRES hdf5_fortran >= 1.14.0"
                 neko_req_libs="$neko_req_libs $HDF5_Fortran_LIBS"
                 neko_req_fcflags="$neko_req_fcflags $HDF5_Fortran_CFLAGS"
              else
                 AC_MSG_NOTICE([HDF5 flags overridden, copied into neko.pc])
              fi
              AC_DEFINE(HAVE_HDF5,[1],
                     [Define if you have the HDF5 library.])
           fi
	fi
])
