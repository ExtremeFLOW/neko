#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "ax_helm_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device HIP Ax
   */
  void hip_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                   void *dxt, void *dyt, void *dzt, void *h1,
                   void *g11, void *g22, void *g33, void *g12,
                   void *g13, void *g23, int *nelv, int *lx) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nelv), 1, 1);

#define CASE(LX)                                                                \
    case LX:                                                                    \
      hipLaunchKernelGGL(HIP_KERNEL_NAME( ax_helm_kernel<real, LX, 1024> ),     \
                         nblcks, nthrds, 0, 0,                                  \
                         (real *) w, (real *) u,                                \
                         (real *) dx, (real *) dy, (real *) dz,                 \
                         (real *) dxt, (real *) dyt, (real *) dzt, (real *) h1, \
                         (real *) g11, (real *) g22, (real *) g33,              \
                         (real *) g12, (real *) g13, (real *) g23);             \
      HIP_CHECK(hipGetLastError());                                             \
      break

    switch(*lx) {
      CASE(2);
      CASE(3);
      CASE(4);
      CASE(5);
      CASE(6);
      CASE(7);
      CASE(8);
      CASE(9);
      CASE(10);
      CASE(11);
      CASE(12);         
    }
  }
}
