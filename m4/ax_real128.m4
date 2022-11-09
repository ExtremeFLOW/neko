#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_REAL128], [

AC_REQUIRE([AC_PROG_FC])
AC_LANG_PUSH([Fortran])

AC_MSG_CHECKING([for REAL128 support])
AC_COMPILE_IFELSE([
       module tt
       use, intrinsic :: iso_fortran_env
       integer, parameter :: qp = REAL128
       end module tt

       program conftest
       use tt       
       real(kind=qp) :: test
       end program conftest
],
[have_real128=yes], [have_real128=no])

AC_LANG_POP([Fortran])

if test "x$have_real128" = "xyes"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_REAL128], [1], [REAL128 support])
else
    AC_MSG_RESULT([no])
fi

])
