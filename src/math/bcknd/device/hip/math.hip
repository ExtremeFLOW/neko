#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include "math_kernel.h"

extern "C" {

  /** Fortran wrapper for copy
   * Copy a vector \f$ a = b \f$
   */
  void hip_copy(void *a, void *b, int *n) {
    hipMemcpyAsync(a, b, (*n) * sizeof(real), hipMemcpyDeviceToDevice);
  }

  /** Fortran wrapper for rzero
   * Zero a real vector
   */
  void hip_rzero(void *a, int *n) {
    hipMemsetAsync(a, 0, (*n) * sizeof(real));
  }

  /** Fortran wrapper for rone
   * Set all elements to one
   */
  void hip_rone(void *a, int *n) {
    hipMemsetAsync(a, 1, (*n) * sizeof(real));
  }

  /** Fortran wrapper for cmult
   * Multiplication by constant c \f$ a = c \cdot a \f$
   */
  void hip_cmult(void *a, real *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cmult_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, *c, *n);
    
  }

  /** Fortran wrapper for cfill
   * Multiplication by constant c \f$ a = c \cdot a \f$
   */
  void hip_cfill(void *a, real *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cfill_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, *c, *n);
    
  }

  /**
   * Fortran wrapper for add2
   * Vector addition \f$ a = a + b \f$
   */
  void hip_add2(void *a, void *b, int *n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, *n);
    
  }
  
  /**
   * Fortran wrapper for add2s1
   * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
   * (multiplication on first argument) 
   */
  void hip_add2s1(void *a, void *b, real *c1, int *n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2s1_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b,
		       *c1, *n);
    
  }

  /**
   * Fortran wrapper for add2s2
   * Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
   * (multiplication on second argument) 
   */
  void hip_add2s2(void *a, void *b, real *c1, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2s2_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b,
		       *c1, *n);

  }

  /**
   * Fortran wrapper for add3s2
   * Vector addition with scalar multiplication \f$ a = c_1 b + c_2 c \f$
   * (multiplication on second argument) 
   */
  void hip_add3s2(void *a, void *b, void *c, real *c1, real *c2, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add3s2_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, (real *) c,
		       *c1, *c2, *n);

  }
  
  /**
   * Fortran wrapper for invcol1
   * Invert a vector \f$ a = 1 / a \f$
   */
  void hip_invcol1(void *a, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(invcol1_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, *n);
  }

  /**
   * Fortran wrapper for invcol2
   * Vector division \f$ a = a / b \f$
   */
  void hip_invcol2(void *a, void *b, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(invcol2_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, *n);
  }
  
  /**
   * Fortran wrapper for col2
   * Vector multiplication with 2 vectors \f$ a = a \cdot b \f$
   */
  void hip_col2(void *a, void *b, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(col2_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, *n);
  }
  
  /**
   * Fortran wrapper for col3
   * Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
   */
  void hip_col3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(col3_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, (real *) c, *n);
  }

  /**
   * Fortran wrapper for subcol3
   * Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
   */
  void hip_subcol3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(subcol3_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, (real *) c, *n);
  }
  
  /**
   * Fortran wrapper for sub2
   * Vector subtraction \f$ a = a - b \f$
   */
  void hip_sub2(void *a, void *b, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sub2_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, *n);
  }
 
  /**
   * Fortran wrapper for sub3
   * Vector subtraction \f$ a = b - c \f$
   */
  void hip_sub3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sub3_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, (real *) c, *n);
  }

  /**
   * Fortran wrapper for addcol3
   * \f$ a = a + b * c \f$
   */
  void hip_addcol3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(addcol3_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, (real *) c, *n);
  }

  /**
   * Fortran wrapper glsc3
   * Weighted inner product \f$ a^T b c \f$
   */
  real hip_glsc3(void *a, void *b, void *c, int *n) {
	
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    real * buf = (real *) malloc(nb * sizeof(real));
    real * buf_d;

    hipMalloc(&buf_d, nb*sizeof(real));
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b,
		       (real *) c, buf_d, *n);

    hipMemcpy(buf, buf_d, nb * sizeof(real), hipMemcpyDeviceToHost);

    real res = 0.0;
    for (int i = 0; i < nb; i++) {
      res += buf[i];
    }

    free(buf);
    hipFree(buf_d);

    return res;
  }

  /**
   * Fortran wrapper glsc2
   * Weighted inner product \f$ a^T b \f$
   */
  real hip_glsc2(void *a, void *b, int *n) {
	
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    real * buf = (real *) malloc(nb * sizeof(real));
    real * buf_d;

    hipMalloc(&buf_d, nb*sizeof(real));
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc2_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b,
		       buf_d, *n);

    hipMemcpy(buf, buf_d, nb * sizeof(real), hipMemcpyDeviceToHost);

    real res = 0.0;
    for (int i = 0; i < nb; i++) {
      res += buf[i];
    }

    free(buf);
    hipFree(buf_d);

    return res;
  }

}
