#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "symmetry_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device symmetry apply vector
   */
  void hip_symmetry_apply_vector(void *xmsk, void *ymsk, void *zmsk,
				 void *x, void *y, void *z,
				 int *m, int *n, int *l) {

    const int max_len = std::max(std::max(*m, *n), *l);
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((max_len) + 1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(symmetry_apply_vector_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (int *) xmsk, (int *) ymsk, (int *) zmsk,
		       (real *) x, (real *) y, (real *) z,
		       *m, *n, *l);
    HIP_CHECK(hipGetLastError());
  }
 
}
