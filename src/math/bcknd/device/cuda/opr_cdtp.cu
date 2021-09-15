#include "cdtp_kernel.h"



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


    switch(*lx) {
    case 2:
      cdtp_kernel<double, 2, 1024>
	<<<nblcks, nthrds>>>((double *) dtx, (double *) x, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dxt, (double *) dyt, (double *) dzt,
			     (double *) B, (double *) jac);
      break;     
    case 3:
      cdtp_kernel<double, 3, 1024>
	<<<nblcks, nthrds>>>((double *) dtx, (double *) x, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dxt, (double *) dyt, (double *) dzt,
			     (double *) B, (double *) jac);
      break;     
    case 4:
      cdtp_kernel<double, 4, 1024>
	<<<nblcks, nthrds>>>((double *) dtx, (double *) x, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dxt, (double *) dyt, (double *) dzt,
			     (double *) B, (double *) jac);
      break;     
    case 5:
      cdtp_kernel<double, 5, 1024>
	<<<nblcks, nthrds>>>((double *) dtx, (double *) x, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dxt, (double *) dyt, (double *) dzt,
			     (double *) B, (double *) jac);
      break;     
    case 6:
      cdtp_kernel<double, 6, 1024>
	<<<nblcks, nthrds>>>((double *) dtx, (double *) x, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dxt, (double *) dyt, (double *) dzt,
			     (double *) B, (double *) jac);
      break;     
    case 7:
      cdtp_kernel<double, 7, 1024>
	<<<nblcks, nthrds>>>((double *) dtx, (double *) x, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dxt, (double *) dyt, (double *) dzt,
			     (double *) B, (double *) jac);
      break;     
    case 8:
      cdtp_kernel<double, 8, 1024>
	<<<nblcks, nthrds>>>((double *) dtx, (double *) x, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dxt, (double *) dyt, (double *) dzt,
			     (double *) B, (double *) jac);
      break;     
    case 9:
      cdtp_kernel<double, 9, 1024>
	<<<nblcks, nthrds>>>((double *) dtx, (double *) x, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dxt, (double *) dyt, (double *) dzt,
			     (double *) B, (double *) jac);
      break;     
    }
  } 
}
