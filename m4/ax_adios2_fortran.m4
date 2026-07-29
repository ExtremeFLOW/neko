AC_DEFUN([AX_ADIOS2_FORTRAN],[
	AC_ARG_WITH([adios2_fortran],
	AS_HELP_STRING([--with-adios2-fortran=DIR],
	[Directory for ADIOS2 with Fortran bindings]),
	[
	if test -d "$withval"; then
	   ac_adios2_fortran_path="$withval";
	fi
	],[with_adios2_fortran=no])

	if test "x${with_adios2_fortran}" != xno; then
	   PATH_SAVED="$PATH"
	   if test -d "$ac_adios2_fortran_path"; then
	      PATH="$ac_adios2_fortran_path/bin:$PATH"
	   fi

	   AC_PATH_PROG(ADIOS2_FORTRAN_CONF,adios2-config,no)

	   if test x"${ADIOS2_FORTRAN_CONF}" = xno; then
	      have_adios2_fortran=no
	      AC_MSG_ERROR([ADIOS2 Fortran support requested but adios2-config was not found])
	   elif ADIOS2_FORTRAN_FCFLAGS=`${ADIOS2_FORTRAN_CONF} --fortran-flags` &&
	        ADIOS2_FORTRAN_LIBS=`${ADIOS2_FORTRAN_CONF} --fortran-libs`; then
	      FCFLAGS="$ADIOS2_FORTRAN_FCFLAGS $FCFLAGS"
	      LIBS="$ADIOS2_FORTRAN_LIBS $LIBS"
              with_adios2_fortran=yes
	      have_adios2_fortran=yes
	      AC_SUBST(have_adios2_fortran)
              AC_DEFINE(HAVE_ADIOS2_FORTRAN,1,[Define if you want to use ADIOS2 with Fortran bindings.])
	   else
	      have_adios2_fortran=no
	      AC_MSG_ERROR([ADIOS2 Fortran support requested but its Fortran configuration is unavailable])
	   fi
            PATH="$PATH_SAVED"

	fi
])
