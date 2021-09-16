#include <hip/hip_runtime.h>

#include "dirichlet_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device dirichlet apply scalar
   */
  void hip_dirichlet_apply_scalar(void *msk, void *x,
				  double *g, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(dirichlet_apply_scalar_kernel<double>),
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (double *) x, *g, *m);
  }
  
  /** 
   * Fortran wrapper for device dirichlet apply vector
   */
  void hip_dirichlet_apply_vector(void *msk, void *x, void *y,
				  void *z, double *g, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(dirichlet_apply_vector_kernel<double>),
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (double *) x,
		       (double *) y, (double *) z, *g, *m);
  }
 
}
