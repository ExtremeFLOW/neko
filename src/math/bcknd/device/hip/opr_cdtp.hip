#include <hip/hip_runtime.h>
#include "cdtp_kernel.h"



extern "C" {

  /** 
   * Fortran wrapper for device hip \f$ D^T X \f$
   */
  void hip_cdtp(void *dtx, void *x,
		void *dr, void *ds, void *dt,
		void *dxt, void *dyt, void *dzt,
		void *B, void *jac, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);


    switch(*lx) {
    case 2:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 2, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;     
    case 3:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 3, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;     
    case 4:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 4, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;     
    case 5:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 5, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;     
    case 6:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 6, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;     
    case 7:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 7, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;     
    case 8:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 8, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;     
    case 9:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 9, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;     
    case 10:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( cdtp_kernel<double, 10, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) dtx, (double *) x, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) B, (double *) jac);
      break;
    }
  } 
}
