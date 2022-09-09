AC_INIT([pyNeko],[0.0.1])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AM_MAINTAINER_MODE
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(Fortran)
AC_PROG_FC(,90)
LT_INIT([disable-static])

PKG_CHECK_MODULES([Neko],[neko >= 0.4.99])
FC=`pkg-config --variable=compiler neko`
FCFLAGS=`pkg-config --cflags neko`
LIBS=`pkg-config --libs neko`
PKG_CHECK_MODULES([JSON_Fortran],[json-fortran >= 8.3.0])

AM_PATH_PYTHON([3.5])
if ! ${HAVE_PYTHON}; then
     AC_MSG_ERROR([Python interpreter not found, pyNeko cannot be installed.])
fi

AC_CONFIG_FILES([Makefile\
                 src/Makefile])

AC_OUTPUT

# Display help text
echo \
"
---------------------------------------------------------
                 _  __ ____ __ __ ____ 
     ___  __ __ / |/ // __// //_// __ \\
    / _ \/ // //    // _/ / ,<  / /_/ /
   / .__/\\_, //_/|_//___//_/|_| \____/ (v${PACKAGE_VERSION})
  /_/   /___/                         

  Python modules will be installed in:

        ${pythondir}

  with \${prefix} = ${prefix} .

  Do not forget to set the PYTHONPATH variable.
---------------------------------------------------------
"
