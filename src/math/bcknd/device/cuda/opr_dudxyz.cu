#include "dudxyz_kernel.h"



extern "C" {

  /** 
   * Fortran wrapper for device cuda derivative kernels
   */
  void cuda_dudxyz(void *du, void *u,
		  void *dr, void *ds, void *dt,
		  void *dx, void *dy, void *dz,
		  void *jacinv, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);


    switch(*lx) {
    case 2:
      dudxyz_kernel<double, 2, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);
      break;     
    case 3:
      dudxyz_kernel<double, 3, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u,
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);

      break;     
    case 4:
      dudxyz_kernel<double, 4, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);

      break;     
    case 5:
      dudxyz_kernel<double, 5, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);

      break;     
    case 6:
      dudxyz_kernel<double, 6, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);

      break;     
    case 7:
      dudxyz_kernel<double, 7, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);

      break;     
    case 8:
      dudxyz_kernel<double, 8, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);

      break;     
    case 9:
      dudxyz_kernel<double, 9, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);

      break;     
    case 10:
      dudxyz_kernel<double, 10, 1024>
	<<<nblcks, nthrds>>>((double *) du, (double *) u, 
			     (double *) dr, (double *) ds, (double *) dt,
			     (double *) dx, (double *) dy, (double *) dz,
			     (double *) jacinv);
      break;
    }
  } 
}
