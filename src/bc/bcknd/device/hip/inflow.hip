#include <hip/hip_runtime.h>

#include "inflow_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device inflow apply vector
   */
  void hip_inflow_apply_vector(void *msk, void *x, void *y,
				  void *z, void *g, int *m) {

    const double gx = ((double *)g)[0];
    const double gy = ((double *)g)[1];
    const double gz = ((double *)g)[2];
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(inflow_apply_vector_kernel<double>),
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (double *) x,
		       (double *) y, (double *) z,
		       gx, gy, gz, *m);
  }
 
}
