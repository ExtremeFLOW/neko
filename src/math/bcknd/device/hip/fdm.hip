#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "fdm_kernel.h"
extern "C" {

  /** Fortran wrapper for tnsr3d **/
  void hip_fdm_do_fast(void *e, void *r, void *s, void *d, int *nl, int *nel) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(*nel, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(fdm_do_fast_kernel<real>),
                       nblcks, nthrds, 0, 0,
                       (real *) e, (real *) r, 
                       (real *) s,(real *) d, *nl);
    HIP_CHECK(hipGetLastError());

  }
}
