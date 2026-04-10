AC_DEFUN([AX_P4EST],[
        AC_ARG_WITH([p4est],
                    AS_HELP_STRING([--with-p4est],
                    [Compile with support for p4est]),
                    [with_p4est=${withval}], [with_p4est=no])

        if test "x${with_p4est}" != xno; then
           PKG_CHECK_MODULES([P4EST],[p4est >= 2.8.7],
                             have_p4est=yes, have_p4est=no)
           if test "x${have_p4est}" = xyes; then
              LIBS="$P4EST_LIBS $LIBS"
              CPPFLAGS="$P4EST_CFLAGS $CPPFLAGS"
              export CPPFLAGS
              AC_DEFINE(HAVE_P4EST,[1],
                        [Define if you have the p4est library.])
           fi
        fi
])
