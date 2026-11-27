AC_PREREQ([2.69])
AC_INIT([gmsh2nmsh],
        m4_esyscmd_s([v=`sed -n 's/^AC_INIT(.neko.,.\(.*\).)$/\1/p' ../../configure.ac 2>/dev/null`; echo "${v:-unknown}"]))
AC_CONFIG_SRCDIR([gmsh2nmsh.F90])
AM_INIT_AUTOMAKE([foreign])

AC_LANG([Fortran])
AC_PROG_FC
AC_FC_PP_SRCEXT([F90])

AC_MSG_CHECKING([whether $FC supports the Fortran 2008 features gmsh2nmsh uses])
AC_LINK_IFELSE([AC_LANG_SOURCE([[
program conftest
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  implicit none
  character(len=:), allocatable :: s
  real(kind=real64) :: x(3)
  integer(kind=int64) :: n
  integer :: u
  s = 'gmsh'
  x = 1.0_real64
  n = size(x, kind=int64)
  open(newunit=u, file='conftest.out', access='stream', &
       form='unformatted', status='replace')
  write(u) int(n, int32), norm2(x), len(s)
  close(u, status='delete')
end program conftest
]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([$FC does not support the Fortran 2008 features gmsh2nmsh needs])])

AC_SUBST([MAN_DATE], [`date +%Y-%m-%d`])

AC_CONFIG_FILES([Makefile gmsh2nmsh.1])
AC_OUTPUT
