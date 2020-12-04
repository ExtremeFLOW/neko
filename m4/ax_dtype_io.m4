#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_DTYPE_IO], [

AC_REQUIRE([AC_PROG_FC])
AC_LANG_PUSH([Fortran])

AC_MSG_CHECKING([for derived type I/O support])
AC_COMPILE_IFELSE([
       module tt
       type :: test_t
       contains
       procedure :: write
       generic :: write(formatted) => write
       end type test_T

       contains

       subroutine write(this, unit, iotype, v_list, iostat, iomsg)
       class(test_t), intent(in) :: this
       integer, intent(in) :: unit
       character(len=*), intent(in) :: iotype
       integer, intent(in) :: v_list(:)
       integer, intent(out) :: iostat
       character(len=*), intent(inout) :: iomsg
       end subroutine write
       end module

       program conftest
       use tt       
       type(test_t) :: t
       end program conftest
],
[ax_dtype_io=yes], [ax_dtype_io=no])

AC_LANG_POP([Fortran])

if test "x$ax_dtype_io" = "xyes"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_DERIVED_TYPE_IO], [1], [Dervied type I/O support])
else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Fortran compiler doesn't support dervied type I/O])
fi

])
