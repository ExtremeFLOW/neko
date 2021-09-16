#include "inflow_kernel.h"
#include <device/device_config.h>

extern "C" {

  /** 
   * Fortran wrapper for device inflow apply vector
   */
  void cuda_inflow_apply_vector(void *msk, void *x, void *y,
				void *z, void *g, int *m) {

    const real gx = ((real *)g)[0];
    const real gy = ((real *)g)[1];
    const real gz = ((real *)g)[2];
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    inflow_apply_vector_kernel<real>
      <<<nblcks, nthrds>>>((int *) msk,
			   (real *) x, (real *) y, (real *) z,
			   gx, gy, gz, *m);
  }
 
}
