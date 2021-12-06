#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include "blasius_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device blasius apply vector
   */
  void hip_blasius_apply_vector(void *msk, void *x, void *y, void *z,
				 void *bla_x, void *bla_y, void *bla_z, int *m) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(blasius_apply_vector_kernel<real>),
		       nblcks, nthrds, 0, 0, (int *) msk, 
			   (real *) x, (real *) y, (real *) z,
			   (real *) bla_x, (real *) bla_y, (real *) bla_z,
			   *m);
  }
 
}
