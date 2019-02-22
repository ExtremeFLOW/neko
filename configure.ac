AC_INIT([neko],[0.0.1])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(Fortran)
AC_PROG_RANLIB
AC_PROG_INSTALL


AC_ARG_ENABLE(contrib,
              AC_HELP_STRING([--enable-contrib],
              [Compile various tools]), 
              [enable_contrib=${enableval}], [enable_contrib=yes])

# Test for a sane fortran environment (^-^)
AC_PROG_FC(,90)
#AC_FC_MODULE_EXTENSION
#AX_F90_MODULE_FLAG
AX_COARRAY
AX_DTYPE_IO


# Checks for Doxygen
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -n "$DOXYGEN"; then
  AC_CHECK_PROGS([DOXYGEN_DEPS], [dot])
fi
AM_CONDITIONAL([ENABLE_DOXYGEN], [test -n "$DOXYGEN_DEPS"])

AC_CONFIG_FILES([Makefile src/Makefile neko.pc])

if test "x${enable_contrib}" = xyes; then
   AC_CONFIG_FILES([contrib/Makefile
		    contrib/rea2nbin/Makefile])
fi
# Doxygen
AC_CONFIG_FILES([doc/Doxyfile doc/Makefile])

AC_OUTPUT
