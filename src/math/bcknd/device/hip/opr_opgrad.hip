#include <hip/hip_runtime.h>
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
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<2, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);
      break;     
    case 3:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<3, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);     
      break;     
    case 4:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<4, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);     
      break;     
    case 5:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<5, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);      
      break;     
    case 6:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<6, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);
      break;     
    case 7:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<7, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);     
      break;     
    case 8:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<8, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);
      break;     
    case 9:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<9, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);
      break;     
    case 10:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( opgrad_kernel<10, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) ux, (double *) uy, (double *) uz, (double *) u,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) drdx, (double *) dsdx, (double *) dtdx,
			 (double *) drdy, (double *) dsdy, (double *) dtdy,
			 (double *) drdz, (double *) dsdz, (double *) dtdz,
			 (double *) w3);
      break;
    }
  } 
}
