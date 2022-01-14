#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "facet_normal_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device facet normal apply surfvec
   */
  void hip_facet_normal_apply_surfvec(void *msk, void *facet,
				      void *x, void *y, void *z,
				      void *u, void *v, void *w,
				      void *nx, void * ny, void *nz,
				      void *area, int *lx, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(facet_normal_apply_surfvec_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (int *) facet,
		       (real *) x, (real *) y, (real *) z,
		       (real *) u, (real *) v, (real *) w,
		       (real *) nx, (real *) ny, (real *) nz,
		       (real *) area, *lx, *m);
    HIP_CHECK(hipGetLastError());
  }
 
}
