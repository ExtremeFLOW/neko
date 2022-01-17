#include "conv1_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>


extern "C" {

  /** 
   * Fortran wrapper for device cuda convective terms
   */
  void cuda_conv1(void *du, void *u,
		  void *vx, void *vy, void *vz,
		  void *dx, void *dy, void *dz,
		  void *drdx, void *dsdx, void *dtdx,
		  void *drdy, void *dsdy, void *dtdy,
		  void *drdz, void *dsdz, void *dtdz,
		  void *jacinv, int *nel, int *gdim, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);

#define CASE(LX)                                                                \
    case LX:                                                                    \
      conv1_kernel<real, LX, 1024>                                              \
	<<<nblcks, nthrds>>>                                                    \
	((real *) du, (real *) u,                                               \
	 (real *) vx, (real *) vy, (real *) vz,                                 \
	 (real *) dx, (real *) dy, (real *) dz,                                 \
	 (real *) drdx, (real *) dsdx, (real *) dtdx,                           \
	 (real *) drdy, (real *) dsdy, (real *) dtdy,                           \
	 (real *) drdz, (real *) dsdz, (real *) dtdz,                           \
	 (real *) jacinv);                                                      \
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
