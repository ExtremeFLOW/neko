AC_PREREQ([2.69])
AC_INIT([neko],[0.7.99])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AM_MAINTAINER_MODE
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AC_LANG(Fortran)
AC_PROG_RANLIB
AC_PROG_INSTALL

AC_ARG_ENABLE(contrib,
              AS_HELP_STRING([--enable-contrib],[Compile various tools]), 
              [enable_contrib=${enableval}], [enable_contrib=yes])

AC_ARG_ENABLE(real,
	      AS_HELP_STRING([--enable-real=Xp  Specify working precision of REAL types:
                    sp -- REAL(kind=4)
	            dp -- REAL(kind=8) (default)
	            qp -- REAL(kind=16)]),
	        [enable_real=${enableval}], [enable_real=dp])

AC_ARG_ENABLE(blk_size,
              AS_HELP_STRING([--enable-blk_size=number of values Specify size of fast memory]),
              [blk_size=${enableval}], [blk_size="10000"])

AC_ARG_ENABLE(device-mpi,
	      AS_HELP_STRING([--enable-device-mpi],[Enable device aware MPI]),
	      [enable_device_mpi=${enableval}], [enable_device_mpi=no])

# Test for a sane fortran environment (^-^)
AC_LANG(Fortran)
AC_PROG_FC(,90)
AX_COARRAY
AX_DTYPE_IO
AX_REAL128

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
   AC_LANG(C++)
   AX_MPI([have_mpi=yes],[have_mpi=no])
   if  test "x${have_mpi}" != xno; then
       CXX="$MPICXX"
   fi
   AC_LANG(Fortran)
else
   AC_MSG_ERROR([Can't find a suitable MPI compiler])   
fi
AX_MPIF08
AX_MPI_DTYPE

AC_FC_PP_SRCEXT([F90])
AX_COMPILER_VENDOR

        
              
# Checks for platform headers
AC_LANG_PUSH([C])
AC_CHECK_HEADERS([cpuid.h])
AC_LANG_POP([C])  

# Set necessary flags to fix issues with various compilers
if test "x$ax_cv_fc_compiler_vendor" = xportland; then
   FCFLAGS="$FCFLAGS  -Mbackslash"
fi

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

# Checks for various platforms
AX_CRAY
AX_HPE_CRAY
AX_SX

# Checks for pFUnit
AX_PFUNIT

# Check for blas and lapack
AS_IF([test "x${is_cray}" = xyes || test "x${is_hpe_cray}" = xyes],
       AX_CRAY_LIBSCI
       if test "x${have_cray_libsci}" = xyes; then
         have_blas_lapack="yes"
       else
         have_blas_lapack="no"
       fi, have_blas_lapack="no")

if test "x${have_blas_lapack}" != xyes; then
  AX_LAPACK
  if test "x${ax_lapack_ok}" != xno; then
     LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
     have_blas_lapack=yes
  else
     have_blas_lapack=no
  fi
fi

if test "x${have_blas_lapack}" != xyes; then
   AC_MSG_ERROR([Can't find a suitable BLAS/LAPACK library])
fi

if test "x${have_darshan}" != xyes; then
    AX_DARSHAN
    if test "x${ax_darshan_ok}" != xno; then
        LIBS="$DARSHAN_LIBS $LIBS"
        have_darshan=yes
    else
        have_darshan=no
    fi
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

# Checks for ADIOS2
AX_ADIOS2

# Checks for gslib
AX_GSLIB

# Select Working precision of REAL
if test "x$enable_real" = xsp; then
   AC_SUBST(NEKO_REAL_TYPE, "sp")
   AC_SUBST(NEKO_C_REAL_TYPE, "c_float")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_REAL")
   AC_SUBST(NEKO_DEV_REAL_TYPE, "float")   				
elif test "x$enable_real" = xdp; then
   AC_SUBST(NEKO_REAL_TYPE, "dp")
   AC_SUBST(NEKO_C_REAL_TYPE, "c_double")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_DOUBLE_PRECISION")
   AC_SUBST(NEKO_DEV_REAL_TYPE, "double")   				
elif test "x$enable_real" = xqp; then
   AC_SUBST(NEKO_REAL_TYPE, "qp")
   AC_SUBST(NEKO_C_REAL_TYPE, "c_long_double")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_REAL16")
   AC_SUBST(NEKO_DEV_REAL_TYPE, "long double")   				
else
   AC_MSG_ERROR([Invalid working precision])
fi

if test "x$have_real128" = xyes; then
   AC_SUBST(NEKO_REAL128_TYPE, "REAL128")
else
   AC_SUBST(NEKO_REAL128_TYPE, "selected_real_kind(1)")
   if test "x$enable_real" = xqp; then
      AC_MSG_ERROR([REAL128 not supported])
   fi
fi                          

AC_SUBST(blk_size, ${blk_size})

# Checks for libxsmm backend
AX_LIBXSMM

# Checks for HIP backend
if test "x${is_hpe_cray}" = xyes; then
   AX_CRAY_ROCM
else
   AX_HIP
fi

# Checks for CUDA backend
if test "x${is_cray}" = xyes; then
  AX_CRAY_CUDATOOLKIT
else
  AX_CUDA
fi

# Checks for OpenCL backend
AX_OPENCL

# Check if device aware MPI is requested
if (test "x${have_cuda}" = xyes ||
    test "x${have_hip}" = xyes  ||
    test "x${have_opencl}" = xyes); then
   if test "x${enable_device_mpi}" = xyes; then

      # Check that necessary accelerator modules are loaded
      if (test "x${is_cray}" = xyes || test "x${is_hpe_cray}" = xyes); then    
      	 AX_CRAY_ACCEL
      fi

      AC_SUBST(device_mpi, .true.)
      AC_DEFINE(HAVE_DEVICE_MPI, 1, [Define if device MPI is requested.])
      
      if test "x${have_cuda}" = xyes; then
         CUDA_CFLAGS="$CUDA_CFLAGS -DHAVE_DEVICE_MPI=1"

         # Make sure we have support for atomicAdd
      	 AS_IF([test "$CUDA_ARCH"],[],[CUDA_ARCH="-arch sm_60"])         
      fi      

      if  test "x${have_hip}" = xyes; then
         HIP_HIPCC_FLAGS="$HIP_HIPCC_FLAGS -DHAVE_DEVICE_MPI=1"
      fi
   else
      AC_SUBST(device_mpi, .false.)
   fi
else
   AC_SUBST(device_mpi, .false.)
fi

# Checks for NCCL backend
AX_NCCL

# Checks for RCCL backend
AX_RCCL

# Checks for NVTX profiling
AX_NVTX

# Checks for ROCTX profiling
AX_ROCTX

# JSON-Fortran
PKG_CHECK_MODULES([JSON_Fortran],[json-fortran >= 7.1.0])
LIBS="$JSON_Fortran_LIBS $LIBS"
FCFLAGS="$JSON_Fortran_CFLAGS $FCFLAGS"

AM_CONDITIONAL([ENABLE_PARMETIS], [test "x$found_domain_decomp" = xyes])
AM_CONDITIONAL([ENABLE_ADIOS2], [test "x${have_adios2}" = xyes])
AM_CONDITIONAL([ENABLE_CUDA], [test "x${have_cuda}" = xyes])
AM_CONDITIONAL([ENABLE_HIP], [test "x${have_hip}" = xyes])
AM_CONDITIONAL([ENABLE_OPENCL], [test "x${have_opencl}" = xyes])

# Set device dependent flags
AC_SUBST(CUDA_ARCH)
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(HIP_HIPCC_FLAGS)

# Set device backend if requested (and found)
if (test "x${have_cuda}" = xyes ||
    test "x${have_hip}" = xyes  ||
    test "x${have_opencl}" = xyes); then
    AC_SUBST(device_bcknd, "1")
else
    AC_SUBST(device_bcknd, "0")
fi

AC_CONFIG_FILES([Makefile\
		 src/Makefile\
		 tests/Makefile\
		 src/comm/comm.F90\
	 	 src/config/neko_config.f90\
		 src/config/num_types.f90\
		 src/device/device_config.h\
		 src/scripts/cltostring.sh\
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
		 tests/space/Makefile\
		 tests/dofmap/Makefile\
		 tests/field/Makefile\
		 tests/mean_field/Makefile\
		 tests/mean_sqr_field/Makefile\
		 tests/gather_scatter/Makefile\
		 tests/krylov/Makefile\
		 tests/bc/Makefile\
		 tests/signal/Makefile\
		 tests/jobctrl/Makefile\
		 tests/device/Makefile\
		 tests/device_math/Makefile\
		 tests/mathops/Makefile\
		 tests/device_mathops/Makefile\
		 tests/octree/Makefile\
                 tests/vector/Makefile\
                 tests/scratch_registry/Makefile\
                 tests/fast3d/Makefile\
                 tests/time_scheme/Makefile\
                 tests/time_scheme_controller/Makefile\
                 tests/time_based_controller/Makefile\
                 tests/point_interpolation/Makefile])

if test "x${enable_contrib}" = xyes; then
   AC_CONFIG_FILES([contrib/Makefile\
		    contrib/rea2nbin/Makefile\
		    contrib/genmeshbox/Makefile\
		    contrib/average_fields_in_time/Makefile\
		    contrib/postprocess_fluid_stats/Makefile\
		    contrib/average_field_in_space/Makefile\
		    contrib/prepart/Makefile])
fi
# Doxygen
AC_CONFIG_FILES([doc/Doxyfile doc/Makefile])

AC_OUTPUT
