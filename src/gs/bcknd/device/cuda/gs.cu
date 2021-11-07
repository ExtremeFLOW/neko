#include <climits>
#include <cstdio>
#include <device/device_config.h>
#include "gs_kernels.h"

#define GS_OP_ADD  1
#define GS_OP_MUL  2
#define GS_OP_MIN  3
#define GS_OP_MAX  4

extern "C" {

  /** 
   * Fortran wrapper for device gather kernels
   */
  void cuda_gather_kernel(void *v, int *m, int *o, void *dg,
			  void *u, int *n, void *gd, int *nb,
			  void *b, void *bo, void *w, int *op) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+ 1024 - 1)/ 1024, 1, 1);

    switch (*op) {
    case GS_OP_ADD:
      cudaMemset(v, 0, (*m) * sizeof(real));
      gather_kernel_add<real>
	<<<nblcks, nthrds>>>((real *) v, *m, *o, (int *) dg,
			     (real *) u, *n, (int *) gd,
			     *nb, (int *) b, (int *) bo);
      break;
    case GS_OP_MUL:
      cudaMemset(v, 1, (*m) * sizeof(real));
      gather_kernel_mul<real>
	<<<nblcks, nthrds>>>((real *) v, *m, *o, (int *) dg,
			     (real *) u, *n, (int *) gd,
			     *nb, (int *) b, (int *) bo);
      break;
    case GS_OP_MIN:
      cudaMemset(v, INT_MAX, (*m) * sizeof(real));
      gather_kernel_min<real>
	<<<nblcks, nthrds>>>((real *) v, *m, *o, (int *) dg,
			     (real *) u, *n, (int *) gd,
			     *nb, (int *) b, (int *) bo);
      break;
    case GS_OP_MAX:
      cudaMemset(v, -INT_MAX, (*m) * sizeof(real));
      gather_kernel_max<real>
	<<<nblcks, nthrds>>>((real *) v, *m, *o, (int *) dg,
			     (real *) u, *n, (int *) gd,
			     *nb, (int *) b, (int *) bo);
      break;
    }
  }

  /**
   * Fortran wrapper for device scatter kernel
   */
  void cuda_scatter_kernel(void *v, int *m, void *dg,
			   void *u, int *n, void *gd,
			   int *nb, void *b, void *bo) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    scatter_kernel<real>
      <<<nblcks, nthrds>>>((real *) v, *m, (int *) dg,
			   (real *) u, *n, (int *) gd,
			   *nb, (int *) b, (int *) bo);
  }
}
