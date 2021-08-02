#include <hip/hip_runtime.h>
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
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<2, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;     
    case 3:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<3, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;     
    case 4:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<4, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;     
    case 5:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<5, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;     
    case 6:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<6, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;     
    case 7:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<7, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;     
    case 8:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<8, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;     
    case 9:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<9, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;     
    case 10:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( dudxyz_kernel<10, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);
      break;
    }
  } 
}
