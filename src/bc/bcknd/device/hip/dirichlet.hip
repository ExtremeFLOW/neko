#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include "dirichlet_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device dirichlet apply scalar
   */
  void hip_dirichlet_apply_scalar(void *msk, void *x,
				  real *g, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(dirichlet_apply_scalar_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (real *) x, *g, *m);
  }
  
  /** 
   * Fortran wrapper for device dirichlet apply vector
   */
  void hip_dirichlet_apply_vector(void *msk, void *x, void *y,
				  void *z, real *g, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(dirichlet_apply_vector_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (real *) x,
		       (real *) y, (real *) z, *g, *m);
  }
 
}
