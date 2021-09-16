#include "cdtp_kernel.h"
#include <device/device_config.h>


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
      cdtp_kernel<real, 2, 1024>
	<<<nblcks, nthrds>>>((real *) dtx, (real *) x, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dxt, (real *) dyt, (real *) dzt,
			     (real *) B, (real *) jac);
      break;     
    case 3:
      cdtp_kernel<real, 3, 1024>
	<<<nblcks, nthrds>>>((real *) dtx, (real *) x, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dxt, (real *) dyt, (real *) dzt,
			     (real *) B, (real *) jac);
      break;     
    case 4:
      cdtp_kernel<real, 4, 1024>
	<<<nblcks, nthrds>>>((real *) dtx, (real *) x, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dxt, (real *) dyt, (real *) dzt,
			     (real *) B, (real *) jac);
      break;     
    case 5:
      cdtp_kernel<real, 5, 1024>
	<<<nblcks, nthrds>>>((real *) dtx, (real *) x, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dxt, (real *) dyt, (real *) dzt,
			     (real *) B, (real *) jac);
      break;     
    case 6:
      cdtp_kernel<real, 6, 1024>
	<<<nblcks, nthrds>>>((real *) dtx, (real *) x, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dxt, (real *) dyt, (real *) dzt,
			     (real *) B, (real *) jac);
      break;     
    case 7:
      cdtp_kernel<real, 7, 1024>
	<<<nblcks, nthrds>>>((real *) dtx, (real *) x, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dxt, (real *) dyt, (real *) dzt,
			     (real *) B, (real *) jac);
      break;     
    case 8:
      cdtp_kernel<real, 8, 1024>
	<<<nblcks, nthrds>>>((real *) dtx, (real *) x, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dxt, (real *) dyt, (real *) dzt,
			     (real *) B, (real *) jac);
      break;     
    case 9:
      cdtp_kernel<real, 9, 1024>
	<<<nblcks, nthrds>>>((real *) dtx, (real *) x, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dxt, (real *) dyt, (real *) dzt,
			     (real *) B, (real *) jac);
      break;     
    }
  } 
}
