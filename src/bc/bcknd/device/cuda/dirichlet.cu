#include "dirichlet_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device dirichlet apply scalar
   */
  void cuda_dirichlet_apply_scalar(void *msk, void *x,
				  double *g, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    dirichlet_apply_scalar_kernel<double>
      <<<nblcks, nthrds>>>((int *) msk, (double *) x, *g, *m);
  }
  
  /** 
   * Fortran wrapper for device dirichlet apply vector
   */
  void cuda_dirichlet_apply_vector(void *msk, void *x, void *y,
				  void *z, double *g, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    dirichlet_apply_vector_kernel<double>
      <<<nblcks, nthrds>>>((int *) msk,
			   (double *) x, (double *) y, (double *) z, *g, *m);
  }
 
}
