#include "dudxyz_kernel.h"
#include <device/device_config.h>


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
      dudxyz_kernel<real, 2, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u,
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);
      break;     
    case 3:
      dudxyz_kernel<real, 3, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u,
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);

      break;     
    case 4:
      dudxyz_kernel<real, 4, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);

      break;     
    case 5:
      dudxyz_kernel<real, 5, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);

      break;     
    case 6:
      dudxyz_kernel<real, 6, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);

      break;     
    case 7:
      dudxyz_kernel<real, 7, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);

      break;     
    case 8:
      dudxyz_kernel<real, 8, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);

      break;     
    case 9:
      dudxyz_kernel<real, 9, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);

      break;     
    case 10:
      dudxyz_kernel<real, 10, 1024>
	<<<nblcks, nthrds>>>((real *) du, (real *) u, 
			     (real *) dr, (real *) ds, (real *) dt,
			     (real *) dx, (real *) dy, (real *) dz,
			     (real *) jacinv);
      break;
    }
  } 
}
