#include "opgrad_kernel.h"



extern "C" {

  /** 
   * Fortran wrapper for device cuda convective terms
   */
  void cuda_opgrad(void *ux, void *uy, void *uz, void *u,
		   void *dx, void *dy, void *dz,
		   void *drdx, void *dsdx, void *dtdx,
		   void *drdy, void *dsdy, void *dtdy,
		   void *drdz, void *dsdz, void *dtdz,
		   void *w3, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);


    switch(*lx) {
    case 2:
      opgrad_kernel<2, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;     
    case 3:
      opgrad_kernel<3, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;     
    case 4:
      opgrad_kernel<4, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;     
    case 5:
      opgrad_kernel<5, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;     
    case 6:
      opgrad_kernel<6, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;     
    case 7:
      opgrad_kernel<7, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;     
    case 8:
      opgrad_kernel<8, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;     
    case 9:
      opgrad_kernel<9, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;     
    case 10:
      opgrad_kernel<10, 1024>
	<<<nblcks, nthrds>>>((double *) ux, (double *) uy, (double *) uz, (double *) u,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) drdx, (double *) dsdx, (double *) dtdx,
			     (double *) drdy, (double *) dsdy, (double *) dtdy,
			     (double *) drdz, (double *) dsdz, (double *) dtdz,
			     (double *) w3);
      break;
    }
  } 
}
