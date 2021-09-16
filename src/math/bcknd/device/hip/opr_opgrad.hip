#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include "opgrad_kernel.h"



extern "C" {

  /** 
   * Fortran wrapper for device hip convective terms
   */
  void hip_opgrad(void *ux, void *uy, void *uz, void *u,
		  void *dx, void *dy, void *dz,
		  void *drdx, void *dsdx, void *dtdx,
		  void *drdy, void *dsdy, void *dtdy,
		  void *drdz, void *dsdz, void *dtdz,
		  void *w3, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);


    switch(*lx) {
    case 2:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 2, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);
      break;     
    case 3:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 3, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);     
      break;     
    case 4:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 4, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);     
      break;     
    case 5:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 5, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);      
      break;     
    case 6:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 6, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);
      break;     
    case 7:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 7, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);     
      break;     
    case 8:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 8, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);
      break;     
    case 9:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 9, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);
      break;     
    case 10:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<real, 10, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (real *) ux, (real *) uy, (real *) uz, (real *) u,
			 (real *) dx, (real *) dy, (real *) dz,
			 (real *) drdx, (real *) dsdx, (real *) dtdx,
			 (real *) drdy, (real *) dsdy, (real *) dtdy,
			 (real *) drdz, (real *) dsdz, (real *) dtdz,
			 (real *) w3);
      break;
    }
  } 
}
