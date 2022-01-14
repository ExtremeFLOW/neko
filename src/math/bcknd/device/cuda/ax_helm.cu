#include "ax_helm_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {

  /** 
   * Fortran wrapper for device CUDA Ax
   */
  void cuda_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                    void *h1, void *g11, void *g22, void *g33, void *g12,
                    void *g13, void *g23, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);

#define CASE(LX)                                                                \
    case LX:                                                                    \
      ax_helm_kernel<real, LX>                                                  \
      <<<nblcks, nthrds>>>((real *) w, (real *) u,                              \
                           (real *) dx, (real *) dy, (real *) dz, (real *) h1,  \
                           (real *) g11, (real *) g22, (real *) g33,            \
                           (real *) g12, (real *) g13, (real *) g23);           \
      CUDA_CHECK(cudaGetLastError());                                           \
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
