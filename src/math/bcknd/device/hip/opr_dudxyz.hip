#include <hip/hip_runtime.h>
#include <device/device_config.h>
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


    switch(*lx) {
    case 2:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 2, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;     
    case 3:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 3, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;     
    case 4:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 4, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;     
    case 5:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 5, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;     
    case 6:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 6, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;     
    case 7:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 7, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;     
    case 8:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 8, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;     
    case 9:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 9, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;     
    case 10:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<real, 10, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) du, (real *) u, 
			 (real *) dr, (real *) ds, (real *) dt,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) jacinv);
      break;
    }
  } 
}
