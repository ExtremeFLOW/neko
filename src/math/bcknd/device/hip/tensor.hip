#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "tensor_kernel.h"
extern "C" {

  /** Fortran wrapper for tnsr3d **/
  void hip_tnsr3d(void *v, int *nv, void *u, int *nu, void *A, void *Bt, void *Ct, int *nel) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(*nel, 1, 1);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(tnsr3d_kernel<real>
),
                       nblcks, nthrds, 0, 0,
                       (real *) v, *nv, (real *) u, 
		       *nu, (real *) A, (real *) Bt, (real *) Ct);
    HIP_CHECK(hipGetLastError());
  }
}
