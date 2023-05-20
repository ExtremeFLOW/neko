#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_CRAY],[
	AC_MSG_CHECKING([for a Cray XT, XE, XC system])
	AC_LANG_PUSH([C])
	AC_EGREP_CPP(yes,
	[#if defined(__CRAYXT) || defined(__CRAYXE) || defined(__CRAYXC)
	  yes
	 #endif
	],
	[AC_MSG_RESULT([yes])
	is_cray="yes"],
	[is_cray="no"
	AC_MSG_RESULT([no])])
	AC_LANG_POP([C])
	AC_SUBST(is_cray)])

AC_DEFUN([AX_HPE_CRAY],[
	AC_MSG_CHECKING([for a HPE Cray system])
	AC_LANG_PUSH([C])
	AC_EGREP_CPP(yes,
	[#if defined(__CRAYXT_COMPUTE_LINUX_TARGET) && (!defined(__CRAYXT) && !defined(__CRAYXE) && !defined(__CRAYXC))
	  yes
	 #endif
	],
	[AC_MSG_RESULT([yes])
	is_hpe_cray="yes"],
	[is_hpe_cray="no"
	AC_MSG_RESULT([no])])
	AC_LANG_POP([C])
        if test x"${is_hpe_cray}" = xyes; then
           AC_MSG_CHECKING([if CrayPE wrappers are used])
           if test "${FC}" = "ftn"; then
              AC_MSG_RESULT([yes])
           else
              is_hpe_cray="no"
              AC_MSG_RESULT([no])
           fi
        fi
      	AC_SUBST(is_hpe_cray)]) 

AC_DEFUN([AX_CRAY_PETSC],[
	AC_MSG_CHECKING([Cray PETSc])
	if test "${CRAY_PETSC_VERSION}"; then
	   have_cray_petsc="yes"
	else
	   have_cray_petsc="no"
	fi
	AC_SUBST(have_cray_petsc)
	if test "x${have_cray_petsc}" = xyes; then
	   AC_DEFINE(HAVE_PETSC,1,[Define if you have the Petsc library.])
	   AC_MSG_RESULT([yes])
	else
	   AC_MSG_RESULT([no])
	fi
])

AC_DEFUN([AX_CRAY_PARMETIS],[
	AC_ARG_WITH([parmetis],[],
		    [
		    ], [with_parmetis=no])
	if test "x${with_parmetis}" != xno; then
	   AC_MSG_CHECKING([Cray ParMETIS])
	   AC_LANG_PUSH([C])
	   AC_EGREP_CPP(yes,
	   [#if defined(__CRAYXC)
	    yes
	   #endif
	   ], [is_cray_xc="yes"], [is_cray_cx="no"])
	   if test "x${is_cray_xc}" = xyes; then
	      if test "${CRAY_TPSL_VERSION}"; then
	      	 have_cray_parmetis="yes"
              else
		have_cray_parmetis="no"
              fi
	   else
	      if test "${CRAY_TRILINOS_VERSION}"; then
	      	 have_cray_parmetis="yes"
	      elif test "${CRAY_PETSC_VERSION}"; then
	      	 have_cray_parmetis="yes"
	      else
	         have_cray_parmetis="no"
	      fi
	   fi
	   AC_SUBST(have_cray_parmetis)
	   if test "x${have_cray_parmetis}" = xyes; then

	     have_parmetis_real64=no
	     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
	     #include <parmetis.h>
	     ],[
	     #if REALTYPEWIDTH != 64
	     #error 'Not 64'
	     #endif
	     ])], [have_parmetis_real64=yes])
	     if test "x${have_parmetis_real64}" = xyes; then
	        AC_DEFINE(HAVE_PARMETIS_REAL64,1, [ParMETIS real_t is 64bit])
	     fi

	     have_parmetis_int64=no
	     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
	     #include <parmetis.h>
	     ],[
	     #if IDXTYPEWIDTH != 64
	     #error 'Not 64'
	     #endif
	     ])], [have_parmetis_int64=yes])
	     if test "x${have_parmetis_int64}" = xyes; then
	        AC_DEFINE(HAVE_PARMETIS_INT64, 1, [ParMETIS idx_t is 64bit])
	     fi

             AC_DEFINE(HAVE_PARMETIS,1,
	         [Define if you have the ParMETIS library])
	     AC_MSG_RESULT([yes])
	  else
	     AC_MSG_RESULT([no])
	  fi
	  AC_LANG_POP([C])
	fi
])

AC_DEFUN([AX_CRAY_ZOLTAN],[
	AC_MSG_CHECKING([Cray Zoltan])
	if test "${CRAY_TRILINOS_VERSION}"; then
	   have_cray_zoltan="yes"
	else
	   have_cray_zoltan="no"
	fi
	AC_SUBST(have_cray_zoltan)
	if test "x${have_cray_zoltan}" = xyes; then
	   AC_DEFINE(HAVE_ZOLTAN,1,
		     [Define if you have the Zoltan library.])
	   AC_MSG_RESULT([yes])
	else
	   AC_MSG_RESULT([no])
	fi
])

AC_DEFUN([AX_CRAY_LIBSCI],[
	AC_MSG_CHECKING([Cray Scientific Libraries])
	if test "${CRAY_LIBSCI_VERSION}"; then
	   AC_MSG_RESULT([yes])
	   have_cray_libsci="yes"
	else
	   AC_MSG_RESULT([no])
	   have_cray_libsci="no"
	fi
	AC_SUBST(have_cray_libsci)
])

AC_DEFUN([AX_CRAY_CUDATOOLKIT],[
	AC_ARG_WITH([cuda],[],
		    [
		    ], [with_cuda=no])
	cuda_bcknd="0"
	if test "x${with_cuda}" != xno; then
	  AC_MSG_CHECKING([Cray CUDA Toolkit])
	  if test "${CRAY_CUDATOOLKIT_VERSION}"; then
	    AC_MSG_RESULT([yes])
	    AC_PATH_PROG(NVCC, nvcc, "no")
	    AS_IF([test "$CUDA_CFLAGS"],[],[CUDA_CFLAGS="-O3"])
	    AC_DEFINE(HAVE_CUDA,1,[Define if you have CUDA.])
	    have_cuda="yes"
	    cuda_bcknd="1"
	  else
	    AC_MSG_RESULT([no])
            AC_MSG_ERROR([Cray CUDA Toolkit not found])
	    have_cuda="no"
	  fi
	fi
	AC_SUBST(cuda_bcknd)
])

AC_DEFUN([AX_CRAY_ROCM],[
	AC_ARG_WITH([hip],[],
		    [
		    ], [with_hip=no])
	hip_bcknd="0"
	if test "x${with_hip}" != xno; then
	  AC_MSG_CHECKING([Cray ROCm Toolkit])
	  if test "${CRAY_ROCM_VERSION}"; then
	    AC_MSG_RESULT([yes])
            AS_IF([test "$HIPCC"],[],[AC_PATH_PROG(HIPCC, hipcc, "no")])
            AS_IF([test "$HIP_HIPCC_FLAGS"],[],[HIP_HIPCC_FLAGS="-O3"])
	    AC_DEFINE(HAVE_HIP,1,[Define if you have HIP.])
	    have_hip="yes"
	    hip_bcknd="1"
	  else
	    AC_MSG_RESULT([no])
            AC_MSG_ERROR([Cray ROCm Toolkit not found])
	    have_hip="no"
	  fi
	fi
	AC_SUBST(hip_bcknd)
])

AC_DEFUN([AX_CRAY_ACCEL], [
	AC_MSG_CHECKING([Cray Accelerator Target])
	if test "${CRAY_ACCEL_TARGET}"; then
	   if test "x${CRAY_ACCEL_TARGET}" = "xhost"; then
	      AC_MSG_RESULT([no])
	      AC_MSG_ERROR([Invalid accelerator target (host)])
	      have_cray_accel="no"	      
	   else
              AC_MSG_RESULT([yes])
	      have_cray_accel="yes"
	   fi
	else
	   AC_MSG_RESULT([no])
	   AC_MSG_ERROR([Cray Accelerator Target not found])
	   have_cray_accel="no"
	fi
	AC_SUBST(have_cray_accel)
])
