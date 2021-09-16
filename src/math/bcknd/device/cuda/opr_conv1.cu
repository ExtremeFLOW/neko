#include "conv1_kernel.h"
#include <device/device_config.h>


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


    switch(*lx) {
    case 2:
      conv1_kernel<real, 2, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;     
    case 3:
      conv1_kernel<real, 3, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;     
    case 4:
      conv1_kernel<real, 4, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;     
    case 5:
      conv1_kernel<real, 5, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;     
    case 6:
      conv1_kernel<real, 6, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;     
    case 7:
      conv1_kernel<real, 7, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;     
    case 8:
      conv1_kernel<real, 8, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;     
    case 9:
      conv1_kernel<real, 9, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;     
    case 10:
      conv1_kernel<real, 10, 1024>
	<<<nblcks, nthrds>>>
	((real *) du, (real *) u,
	 (real *) vx, (real *) vy, (real *) vz,
	 (real *) dx, (real *) dy, (real *) dz,
	 (real *) drdx, (real *) dsdx, (real *) dtdx,
	 (real *) drdy, (real *) dsdy, (real *) dtdy,
	 (real *) drdz, (real *) dsdz, (real *) dtdz,
	 (real *) jacinv);
      break;
    }
  } 
}
