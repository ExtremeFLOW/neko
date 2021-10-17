#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_MPI_DTYPE], [

AC_REQUIRE([AC_PROG_FC])
AC_LANG_PUSH([Fortran])

AC_MSG_CHECKING([if MPI datatypes can be parameters])
AC_COMPILE_IFELSE([
       program conftest
       use mpi_f08
       type(MPI_Datatype), parameter :: MPI_DTYPE = MPI_REAL
       end program conftest],
[ax_mpi_const_dtype=yes], [ax_mpi_const_dtype=no])

AC_LANG_POP([Fortran])

if test "x$ax_mpi_const_dtype" = "xyes"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_MPI_PARAM_DTYPE], [1], [MPI supports datatypes as parameters])
else
    AC_MSG_RESULT([no])
fi

])


