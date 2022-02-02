#include <device/device_config.h>
#include <device/cuda/check.h>
#include "tensor_kernel.h"
extern "C" {

  /** Fortran wrapper for tnsr3d **/
  void cuda_tnsr3d(void *v, int *nv, void *u, int *nu, void *A, void *Bt, void *Ct, int *nel) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(*nel, 1, 1);
    
    tnsr3d_kernel<real>
      <<<nblcks, nthrds>>>((real *) v, *nv, (real *) u, 
			   *nu, (real *) A, (real *) Bt, (real *) Ct);
    CUDA_CHECK(cudaGetLastError());
    
  }
}
