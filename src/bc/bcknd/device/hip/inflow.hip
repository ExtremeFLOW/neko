#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include "inflow_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device inflow apply vector
   */
  void hip_inflow_apply_vector(void *msk, void *x, void *y,
				  void *z, void *g, int *m) {

    const real gx = ((real *)g)[0];
    const real gy = ((real *)g)[1];
    const real gz = ((real *)g)[2];
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(inflow_apply_vector_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (real *) x,
		       (real *) y, (real *) z,
		       gx, gy, gz, *m);
  }
 
}
