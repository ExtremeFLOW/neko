#include <hip/hip_runtime.h>

#include "no_slip_wall_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device no-slop wall apply vector
   */
  void hip_no_slip_wall_apply_scalar(void *msk, void *x, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(no_slip_wall_apply_scalar_kernel<double>),
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (double *) x, *m);
  }
  
  /** 
   * Fortran wrapper for device no-slop wall apply vector
   */
  void hip_no_slip_wall_apply_vector(void *msk, void *x, void *y,
				     void *z, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(no_slip_wall_apply_vector_kernel<double>),
		       nblcks, nthrds, 0, 0, (int *) msk,
		       (double *) x, (double *) y, (double *) z, *m); 
  }
 
}
