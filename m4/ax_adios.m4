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

	   AC_CHECK_PROG(ADIOS2CONF,adios2-config,yes)

	   if test x"${ADIOS2CONF}" == x"yes"; then
	      ADIOS2_CXXFLAGS=`adios2-config --cxx-flags`
	      CXXFLAGS="$ADIOS2_CXXFLAGS $CXXFLAGS"

	      ADIOS2_LDFLAGS=`adios2-config --cxx-libs`
	      LDFLAGS="$ADIOS2_LDFLAGS $LDFLAGS"
      	      with_adios2=yes
	      have_adios2=yes
	      AC_SUBST(have_adios2)
	    else
	      with_adios2=no
	    fi
            PATH="$PATH_SAVED"
	    
	fi
])
	
