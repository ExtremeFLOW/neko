AC_DEFUN([AX_ADIOS2],[
	AC_ARG_WITH([adios2],
	AS_HELP_STRING([--with-adios2=DIR],
	[Directory for ADIOS2]),
	[
	if test -d "$withval"; then
	   ac_adios2_path="$withval";
	fi		
	],[with_adios2=no])

	if test "x${with_adios2}" != xno; then
	   PATH_SAVED="$PATH"
	   if test -d "$ac_adios2_path"; then
	      PATH="$ac_adios2_path/bin:$PATH"
	   fi

	   AC_PATH_PROG(ADIOS2CONF,adios2-config,no)

	   if test x"${ADIOS2CONF}" = xno; then
	      have_adios2=no
	      AC_MSG_ERROR([ADIOS2 requested but adios2-config was not found])
	   elif ADIOS2_CXXFLAGS=`${ADIOS2CONF} --cxx-flags` &&
	        ADIOS2_LDFLAGS=`${ADIOS2CONF} --cxx-libs`; then
	      CXXFLAGS="$ADIOS2_CXXFLAGS $CXXFLAGS"
	      LDFLAGS="$ADIOS2_LDFLAGS $LDFLAGS"
      	      with_adios2=yes
	      have_adios2=yes
	      AC_SUBST(have_adios2)
              AC_DEFINE(HAVE_ADIOS2,1,[Define if you have ADIOS2.])
	   else
	      have_adios2=no
	      AC_MSG_ERROR([ADIOS2 requested but its C++ configuration is unavailable])
	   fi
            PATH="$PATH_SAVED"
	    
	fi
])
	
