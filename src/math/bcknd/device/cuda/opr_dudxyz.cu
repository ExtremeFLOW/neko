#include "dudxyz_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {

  /** 
   * Fortran wrapper for device cuda derivative kernels
   */
  void cuda_dudxyz(void *du, void *u,
                  void *dr, void *ds, void *dt,
                  void *dx, void *dy, void *dz,
                  void *jacinv, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);

#define CASE(LX)                                                                \
    case LX:                                                                    \
      dudxyz_kernel<real, 2, 1024>                                              \
        <<<nblcks, nthrds>>>((real *) du, (real *) u,                           \
                             (real *) dr, (real *) ds, (real *) dt,             \
                             (real *) dx, (real *) dy, (real *) dz,             \
                             (real *) jacinv);                                  \
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
    }
  } 
}
