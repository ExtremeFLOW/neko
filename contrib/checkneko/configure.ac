AC_PREREQ([2.69])
AC_INIT([checkneko],
        m4_esyscmd_s([v=`sed -n 's/^AC_INIT(.neko.,.\(.*\).)$/\1/p' ../../configure.ac 2>/dev/null`; echo "${v:-unknown}"]))
AC_CONFIG_SRCDIR([checkneko.F90])
AC_CONFIG_MACRO_DIR([../../m4])
AM_INIT_AUTOMAKE([foreign])

AC_LANG([Fortran])
AC_PROG_FC
AC_FC_PP_SRCEXT([F90])

# checkneko needs json-fortran, which Neko's configure has already found.
# Its modules must have been built by the same compiler as FC, otherwise the
# tool is skipped rather than failing the build of Neko.
PKG_PROG_PKG_CONFIG
PKG_CHECK_MODULES([JSON_Fortran], [json-fortran >= 8.3.0],
                  [have_json_fortran=yes], [have_json_fortran=no])
if test "x${have_json_fortran}" = xyes; then
   AC_MSG_CHECKING([whether $FC can use json-fortran])
   save_FCFLAGS="$FCFLAGS"
   save_LIBS="$LIBS"
   FCFLAGS="$FCFLAGS $JSON_Fortran_CFLAGS"
   LIBS="$JSON_Fortran_LIBS $LIBS"
   AC_LINK_IFELSE([AC_LANG_SOURCE([[
program conftest
  use json_module, only : json_core, json_value
  implicit none
  type(json_core) :: core
  type(json_value), pointer :: p
  call core%initialize()
  call core%deserialize(p, '{"case": {}}')
  call core%destroy(p)
end program conftest
]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      have_json_fortran=no])
   FCFLAGS="$save_FCFLAGS"
   LIBS="$save_LIBS"
fi
if test "x${have_json_fortran}" != xyes; then
   AC_MSG_WARN([json-fortran cannot be used with $FC, checkneko is not built])
fi
AM_CONDITIONAL([BUILD_CHECKNEKO], [test "x${have_json_fortran}" = xyes])

AC_SUBST([MAN_DATE], [`date +%Y-%m-%d`])

AC_CONFIG_FILES([Makefile checkneko.1])
AC_OUTPUT
