#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_COARRAY], [

AC_REQUIRE([AC_PROG_FC])
AC_LANG_PUSH([Fortran])

AC_MSG_CHECKING([for coarray support])
AC_COMPILE_IFELSE([
       program conftest
       integer :: i[[*]]
       integer :: image_num
       image_num = this_image()
       sync all
       sync images(*)
       end program conftest],
[ax_coarray=yes], [ax_coarray=no])

AC_LANG_POP([Fortran])

if test "x$ax_coarray" = "xyes"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_COARRAY], [1], [Coarray Fortran support])
else
    AC_MSG_RESULT([no])
fi

])
