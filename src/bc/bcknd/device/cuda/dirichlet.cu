#include "dirichlet_kernel.h"
#include <device/device_config.h>

extern "C" {

  /** 
   * Fortran wrapper for device dirichlet apply scalar
   */
  void cuda_dirichlet_apply_scalar(void *msk, void *x,
				  real *g, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    dirichlet_apply_scalar_kernel<real>
      <<<nblcks, nthrds>>>((int *) msk, (real *) x, *g, *m);
  }
  
  /** 
   * Fortran wrapper for device dirichlet apply vector
   */
  void cuda_dirichlet_apply_vector(void *msk, void *x, void *y,
				  void *z, real *g, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    dirichlet_apply_vector_kernel<real>
      <<<nblcks, nthrds>>>((int *) msk,
			   (real *) x, (real *) y, (real *) z, *g, *m);
  }
 
}
