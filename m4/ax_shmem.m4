#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_SHMEM],[
	AC_ARG_WITH([openshmem],
		    AS_HELP_STRING([--with-openshmem=DIR],
		    [Compile with support for OpenSHMEM]),
		    [
		    if test -d "$withval"; then
		       ac_shmem_path="$withval";
		       AS_IF([test -d "$ac_shmem_path/lib64"],
			     [suffix="64"],[suffix=""])
		       SHMEM_LDFLAGS="-L$ac_shmem_path/lib$suffix"
		    fi
		    ], [with_openshmem=no])
	have_openshmem=no
	if test "x${with_openshmem}" != xno; then
	   LDFLAGS_SAVED="$LDFLAGS"
	   LDFLAGS="$SHMEM_LDFLAGS $LDFLAGS"
	   export LDFLAGS

	   AC_LANG_PUSH([C])

	   # Sandia OpenSHMEM (SOS) provides libsma, whereas the OSHMEM
	   # component of Open MPI provides liboshmem. Neko only uses the
	   # C API (through iso_c_binding), so no Fortran bindings are
	   # needed from the implementation.
	   AC_CHECK_LIB([sma], [shmem_init],
			[have_openshmem=yes; SHMEM_LIBS="-lsma"],
			[AC_CHECK_LIB([oshmem], [shmem_init],
				      [have_openshmem=yes
				       SHMEM_LIBS="-loshmem"],
				      [have_openshmem=no], [])], [])

	   # The gather-scatter backend relies on the signaling operations
	   # (put with signal) added in OpenSHMEM 1.5, which are missing in
	   # implementations that only provide 1.4
	   if test "x${have_openshmem}" = xyes; then
	      LIBS_SAVED="$LIBS"
	      LIBS="$SHMEM_LIBS $LIBS"
	      have_shmem_signal=yes
	      AC_CHECK_FUNC([shmem_putmem_signal_nbi], [],
			    [have_shmem_signal=no])
	      AC_CHECK_FUNC([shmem_signal_wait_until], [],
			    [have_shmem_signal=no])
	      LIBS="$LIBS_SAVED"
	   fi

	   AC_LANG_POP([C])

	   if test "x${have_openshmem}" != xyes; then
	      AC_MSG_ERROR([OpenSHMEM requested but neither libsma (SOS) nor liboshmem (Open MPI) was found])
	   elif test "x${have_shmem_signal}" != xyes; then
	      AC_MSG_ERROR([OpenSHMEM found, but it lacks the signaling operations (OpenSHMEM 1.5) required by Neko])
	   else
	      AC_DEFINE(HAVE_OPENSHMEM, 1, [Define if you have OPENSHMEM.])
	      LIBS="$SHMEM_LIBS $LIBS"
	      LDFLAGS="$SHMEM_LDFLAGS $LDFLAGS_SAVED"
	   fi
	fi
	AC_SUBST(have_openshmem)
])
