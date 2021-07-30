AC_INIT([neko],[0.2.x])
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

AC_ARG_ENABLE(real,
	      AC_HELP_STRING([--enable-real=Xp  Specify working precision of REAL types:
                    sp -- REAL(kind=4)
	            dp -- REAL(kind=8) (default)
	            qp -- REAL(kind=16)]),
	        [enable_real=${enableval}], [enable_real=dp])

AC_ARG_ENABLE(hip,
	      AC_HELP_STRING([--enable-hip],[Compile with HIP backend]),
	                     [enable_hip=${enableval}], [enable_hip=no])

AC_ARG_ENABLE(cuda,
	      AC_HELP_STRING([--enable-cuda],[Compile with CUDA backend]),
	                     [enable_cuda=${enableval}], [enable_cuda=no])

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
   AC_LANG(C)
   AX_MPI([have_mpi=yes],[have_mpi=no])
   if  test "x${have_mpi}" != xno; then
       CC="$MPICC"
   else
       AC_MSG_ERROR([Can't find a suitable C MPI compiler])
   fi
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

# Checks for flinter
AC_CHECK_PROGS([FLINT], [flint])
AM_CONDITIONAL([ENABLE_FLINT], [test -n "$FLINT"])

AX_CRAY
AX_SX

# Checks for pFUnit
AX_PFUNIT

# Check for blas and lapack
AX_LAPACK
if test "x${ax_lapack_ok}" != xno; then
   LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
else
   AC_MSG_ERROR([Can't find a suitable BLAS/LAPACK library])
fi

found_domain_decomp="no"

# Checks for ParMetis
if test "x${is_cray}" = xyes; then
   AX_CRAY_PARMETIS
   if test "x${have_cray_parmetis}" = xyes; then
      found_domain_decomp="yes"
   fi
fi

if test "x${found_domain_decomp}" != xyes; then
   AX_PARMETIS
   if test "${PARMETIS_LIBS}"; then
       LIBS="$PARMETIS_LIBS $LIBS"
       found_domain_decomp="yes"
   fi
fi

# Select Working precision of REAL
if test "x$enable_real" = xsp; then
   AC_SUBST(NEKO_REAL_TYPE, "sp")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_REAL")
elif test "x$enable_real" = xdp; then
   AC_SUBST(NEKO_REAL_TYPE, "dp")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_DOUBLE_PRECISION")
elif test "x$enable_real" = xqp; then
   AC_SUBST(NEKO_REAL_TYPE, "qp")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_REAL16")
else
   AC_MSG_ERROR([Invalid working precision])
fi

# Checks for libxsmm backend
AX_LIBXSMM

# Checks for HIP backend
AX_HIP

# Checks for CUDA backend
AX_CUDA


AM_CONDITIONAL([ENABLE_PARMETIS], [test "x$found_domain_decomp" = xyes])
AM_CONDITIONAL([ENABLE_CUDA], [test "x${have_cuda}" = xyes])
AM_CONDITIONAL([ENABLE_HIP], [test "x${have_hip}" = xyes])



AC_CONFIG_FILES([Makefile\
		 src/Makefile\
		 tests/Makefile\
		 src/comm.f90\
	 	 src/neko_config.f90\
		 src/num_types.f90\
		 makeneko\
		 neko.pc])

# Config tests
AC_CONFIG_FILES([tests/stack/Makefile\
		 tests/tuple/Makefile\
 		 tests/point/Makefile\
		 tests/htable/Makefile\
		 tests/uset/Makefile\
		 tests/quad/Makefile\
		 tests/hex/Makefile\
		 tests/tet/Makefile\
		 tests/tri/Makefile\
		 tests/math/Makefile\
		 tests/mesh/Makefile\
		 tests/dofmap/Makefile\
		 tests/gather_scatter/Makefile\
		 tests/krylov/Makefile\
		 tests/bc/Makefile\
		 tests/signal/Makefile\
		 tests/jobctrl/Makefile])

if test "x${enable_contrib}" = xyes; then
   AC_CONFIG_FILES([contrib/Makefile\
		    contrib/rea2nbin/Makefile\
		    contrib/prepart/Makefile])
fi
# Doxygen
AC_CONFIG_FILES([doc/Doxyfile doc/Makefile])

AC_OUTPUT
