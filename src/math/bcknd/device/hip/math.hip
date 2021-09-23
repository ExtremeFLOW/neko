#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include "math_kernel.h"

extern "C" {

  /** Fortran wrapper for copy
   * Copy a vector \f$ a = b \f$
   */
  void hip_copy(void *a, void *b, int *n) {
    hipMemcpy(a, b, (*n) * sizeof(real), hipMemcpyDeviceToDevice);
  }

  /** Fortran wrapper for rzero
   * Zero a real vector
   */
  void hip_rzero(void *a, int *n) {
    hipMemset(a, 0, (*n) * sizeof(real));
  }

  
  /**
   * Fortran wrapper for add2s1
   * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
   * (multiplication on first argument) 
   */
  void hip_add2s1(void *a, void *b, real *c1, int *n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(add2s1_kernel,
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

    hipLaunchKernelGGL(add2s2_kernel,
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b,
		       *c1, *n);

  }

  /**
   * Fortran wrapper for invcol2
   * Vector division \f$ a = a / b \f$
   */
  void hip_invcol2(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(invcol2_kernel,
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

    hipLaunchKernelGGL(col2_kernel,
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

    hipLaunchKernelGGL(col3_kernel,
		       nblcks, nthrds, 0, 0,
		       (real *) a, (real *) b, (real *) c, *n);
  }
  

  /**
   * Fortran wrapper for sub3
   * Vector subtraction \f$ a = b - c \f$
   */
  void hip_sub3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(sub3_kernel,
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

    hipLaunchKernelGGL(addcol3_kernel,
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
     
    hipLaunchKernelGGL(glsc3_kernel,
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
     
    hipLaunchKernelGGL(glsc2_kernel,
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
