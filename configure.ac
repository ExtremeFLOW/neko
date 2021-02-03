AC_INIT([neko],[0.0.1])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AC_LANG(Fortran)
AC_PROG_RANLIB
AC_PROG_INSTALL

AC_ARG_ENABLE(contrib,
              AC_HELP_STRING([--enable-contrib],
              [Compile various tools]), 
              [enable_contrib=${enableval}], [enable_contrib=yes])

# Test for a sane fortran environment (^-^)
AC_LANG(Fortran)
AC_PROG_FC(,90)
AX_COARRAY
AX_DTYPE_IO

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

AC_FC_PP_SRCEXT([F90])
AX_COMPILER_VENDOR

# Store build information (date, host, FC)
AC_SUBST(NEKO_BUILD_INFO,
        ["(build: $(date +%Y-%m-%d) on $host using $ax_cv_fc_compiler_vendor)"], [])

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

# Check for blas and lapack
AX_LAPACK
if test "x${ax_lapack_ok}" != xno; then
   LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
else
   AC_MSG_ERROR([Can't find a suitable BLAS/LAPACK library])
fi


AC_CONFIG_FILES([Makefile\
		 src/Makefile\
		 tests/Makefile\
	 	 src/neko_config.f90\
		 makenek\
		 neko.pc])

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
		 tests/gather_scatter/Makefile\
		 tests/krylov/Makefile\
		 tests/bc/Makefile])

if test "x${enable_contrib}" = xyes; then
   AC_CONFIG_FILES([contrib/Makefile
		    contrib/rea2nbin/Makefile])
fi
# Doxygen
AC_CONFIG_FILES([doc/Doxyfile doc/Makefile])

AC_OUTPUT
