#include "conv1_kernel.h"



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
      conv1_kernel<2, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;     
    case 3:
      conv1_kernel<3, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;     
    case 4:
      conv1_kernel<4, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;     
    case 5:
      conv1_kernel<5, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;     
    case 6:
      conv1_kernel<6, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;     
    case 7:
      conv1_kernel<7, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;     
    case 8:
      conv1_kernel<8, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;     
    case 9:
      conv1_kernel<9, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;     
    case 10:
      conv1_kernel<10, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) vx, (double *) vy, (double *) vz,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) jacinv);
      break;
    }
  } 
}
