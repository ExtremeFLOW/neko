#include "cdtp_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>


extern "C" {

  /** 
   * Fortran wrapper for device cuda \f$ D^T X \f$
   */
  void cuda_cdtp(void *dtx, void *x,
                 void *dr, void *ds, void *dt,
                 void *dxt, void *dyt, void *dzt,
                 void *B, void *jac, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);

#define CASE(LX)                                                                \
    case LX:                                                                    \
      cdtp_kernel<real, LX, 1024>                                               \
        <<<nblcks, nthrds>>>((real *) dtx, (real *) x,                          \
                             (real *) dr, (real *) ds, (real *) dt,             \
                             (real *) dxt, (real *) dyt, (real *) dzt,          \
                             (real *) B, (real *) jac);                         \
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
