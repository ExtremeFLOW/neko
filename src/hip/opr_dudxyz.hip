#include <hip/hip_runtime.h>


__global__ void dudxyz_kernel10(double * __restrict__ du,
				const double * __restrict__ u,
				const double * __restrict__ dr,
				const double * __restrict__ ds,
				const double * __restrict__ dt,
				const double * __restrict__ dx,
				const double * __restrict__ dy,
				const double * __restrict__ dz,
				const double * __restrict__ jacinv);


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
    case 10:
      hipLaunchKernelGGL(dudxyz_kernel10,
			 nblcks, nthrds, 0, 0,
			 (double *) du, (double *) u, 
			 (double *) dr, (double *) ds, (double *) dt,
			 (double *) dx, (double *) dy, (double *) dz,
			 (double *) jacinv);      
    }
  } 
}
