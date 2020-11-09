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

AX_BLAS

# Test for a working MPI compiler
AX_MPI([have_mpi=yes], [have_mpi=no])
if test "x${have_mpi}" != xno; then
   FC="$MPIFC"
   LIBS="$MPILIBS $LIBS"
#   AC_LANG(C)
#   AX_MPI([have_mpi=yes],[have_mpi=no])
#   if  test "x${have_mpi}" != xno; then
#       CC="$MPICC"
#   else
#       AC_MSG_ERROR([Can't find a suitable C MPI compiler])
#   fi
   AC_LANG(Fortran)
else
   AC_MSG_ERROR([Can't find a suitable MPI compiler])   
fi

#AX_PARMETIS
#if test "${PARMETIS_LIBS}"; then
#   LIBS="$PARMETIS_LIBS $LIBS"
#fi			

#AX_GSLIB
#if test "${GSLIB_LIBS}"; then
#   LIBS="$GSLIB_LIBS $LIBS"
#else
#   AC_MSG_ERROR([Can't find a suitable Gather-Scatter library])
#fi			


# Checks for Doxygen
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -n "$DOXYGEN"; then
  AC_CHECK_PROGS([DOXYGEN_DEPS], [dot])
fi
AM_CONDITIONAL([ENABLE_DOXYGEN], [test -n "$DOXYGEN_DEPS"])

# Checks for makedepf90
AC_CHECK_PROGS([MAKEDEPF90], [makedepf90])
AM_CONDITIONAL([ENABLE_MAKEDEPF90], [test -n "$MAKEDEPF90"])

# Checks for pFUnit
AX_PFUNIT


AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile neko.pc])

# Config tests
AC_CONFIG_FILES([tests/stack/Makefile\
		 tests/tuple/Makefile\
 		 tests/point/Makefile\
		 tests/htable/Makefile\
		 tests/uset/Makefile\
		 tests/quad/Makefile\
		 tests/hex/Makefile\
		 tests/math/Makefile\
		 tests/mesh/Makefile\
		 tests/dofmap/Makefile\
		 tests/krylov/Makefile])

if test "x${enable_contrib}" = xyes; then
   AC_CONFIG_FILES([contrib/Makefile
		    contrib/rea2nbin/Makefile])
fi
# Doxygen
AC_CONFIG_FILES([doc/Doxyfile doc/Makefile])

AC_OUTPUT
