#include <hip/hip_runtime.h>
#include "pc_jacobi_kernel.h"



extern "C" {

  /** 
   * Fortran wrapper for device hip pc jacobi update
   */
  void hip_jacobi_update(void *d, 
			 void *dxt, void *dyt, void *dzt,
			 void *G11, void *G22, void *G33,
			 void *G12, void *G13, void *G23,
			 int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);


    switch(*lx) {
    case 2:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<2, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;     
    case 3:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<3, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;     
    case 4:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<4, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;     
    case 5:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<5, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;     
    case 6:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<6, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;     
    case 7:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<7, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;     
    case 8:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<8, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;     
    case 9:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<9, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;     
    case 10:
      hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<10, 1024> ),
			 nblcks, nthrds, 0, 0,
			 (double *) d, 
			 (double *) dxt, (double *) dyt, (double *) dzt,
			 (double *) G11, (double *) G22, (double *) G33,
			 (double *) G12, (double *) G13, (double *) G23);
      break;
    }
  } 
}
