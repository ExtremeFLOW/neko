#include "no_slip_wall_kernel.h"
#include <device/device_config.h>

extern "C" {

  /** 
   * Fortran wrapper for device no-slop wall apply vector
   */
  void cuda_no_slip_wall_apply_scalar(void *msk, void *x, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    no_slip_wall_apply_scalar_kernel<real>
      <<<nblcks, nthrds>>>((int *) msk, (real *) x, *m);
  }
  
  /** 
   * Fortran wrapper for device no-slop wall apply vector
   */
  void cuda_no_slip_wall_apply_vector(void *msk, void *x, void *y,
				     void *z, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    no_slip_wall_apply_vector_kernel<real>
      <<<nblcks, nthrds>>>((int *) msk,
			   (real *) x, (real *) y, (real *) z, *m);
  }
 
}
