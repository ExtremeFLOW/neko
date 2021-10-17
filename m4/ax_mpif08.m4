#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_MPIF08], [

AC_REQUIRE([AC_PROG_FC])
AC_LANG_PUSH([Fortran])

AC_MSG_CHECKING([for Fortran 2008 MPI support])
AC_COMPILE_IFELSE([
       program conftest
       use mpi_f08
       end program conftest],
[ax_mpif08=yes], [ax_mpif08=no])

AC_LANG_POP([Fortran])

if test "x$ax_mpif08" = "xyes"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_MPI_F08], [1], [Fortran 2008 MPI support])
else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([MPI compiler doesn't have a Fortran 2008 module])
fi

])


