#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "dudxyz_kernel.h"



extern "C" {

  /** 
   * Fortran wrapper for device hip derivative kernels
   */
  void hip_dudxyz(void *du, void *u,
		  void *dr, void *ds, void *dt,
		  void *dx, void *dy, void *dz,
		  void *jacinv, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);

#define CASE(LX)                                                                \
    case LX:                                                                    \
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, LX, 1024> ),      \
			 nblcks, nthrds, 0, 0,                                  \
			 (real *) du, (real *) u,                               \
			 (real *) dr, (real *) ds, (real *) dt,                 \
			 (real *) dx, (real *) dy, (real *) dz,                 \
			 (real *) jacinv);                                      \
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
    }
  } 
}
