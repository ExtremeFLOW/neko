#include "ax_helm_kernel.h"
#include <device/device_config.h>

extern "C" {

  /** 
   * Fortran wrapper for device CUDA Ax
   */
  void cuda_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
		    void *h1, void *g11, void *g22, void *g33, void *g12,
		    void *g13, void *g23, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);

    switch(*lx) {
    case 2:
      {
	ax_helm_kernel<real, 2>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 3:
      {
	ax_helm_kernel<real, 3>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 4:
      {
	ax_helm_kernel<real, 4>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 5:
      {
	ax_helm_kernel<real, 5>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 6:
      {
	ax_helm_kernel<real, 6>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 7:
      {
	ax_helm_kernel<real, 7>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 8:
      {
	ax_helm_kernel<real, 8>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 9:
      {
	ax_helm_kernel<real, 9>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 10:
      {
	ax_helm_kernel<real, 10>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 11:
      {
	ax_helm_kernel<real, 11>
	<<<nblcks, nthrds>>>((real *) w, (real *) u,
			     (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			     (real *) g11, (real *) g22, (real *) g33,
			     (real *) g12, (real *) g13, (real *) g23);
      }
      break;
    case 12:
      {
	ax_helm_kernel<real, 12>
	  <<<nblcks, nthrds>>>((real *) w, (real *) u,
			       (real *) dx, (real *) dy, (real *) dz, (real *) h1,
			       (real *) g11, (real *) g22, (real *) g33,
			       (real *) g12, (real *) g13, (real *) g23);
      }
    }    
  }
}
