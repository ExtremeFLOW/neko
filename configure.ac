AC_INIT([neko],[0.0.1])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AC_LANG(Fortran)
AC_CONFIG_MACRO_DIR([m4])

# Test for a sane fortran environment (^-^)
AC_PROG_FC(,90)
#AC_FC_MODULE_EXTENSION
#AX_F90_MODULE_FLAG
AX_COARRAY

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
