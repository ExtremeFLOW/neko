#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "cdtp_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device hip \f$ D^T X \f$
   */
  void hip_cdtp(void *dtx, void *x,
                void *dr, void *ds, void *dt,
                void *dxt, void *dyt, void *dzt,
                void *B, void *jac, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);

#define CASE(LX)                                                                \
    case LX:                                                                    \
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<real, LX, 1024> ),        \
                         nblcks, nthrds, 0, 0,                                  \
                         (real *) dtx, (real *) x,                              \
                         (real *) dr, (real *) ds, (real *) dt,                 \
                         (real *) dxt, (real *) dyt, (real *) dzt,              \
                         (real *) B, (real *) jac);                             \
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
