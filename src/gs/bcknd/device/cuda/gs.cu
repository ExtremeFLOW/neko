#include <climits>
#include <cstdio>

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
			  void *u, int *n, void *gd, void *w, int *op) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+ 1024 - 1)/ 1024, 1, 1);

    switch (*op) {
    case GS_OP_ADD:
      cudaMemset(v, 0, (*m) * sizeof(double));
      gather_kernel_add<double>
	<<<nblcks, nthrds>>>((double *) v, *m, *o, (int *) dg,
			     (double *) u, *n, (int *) gd, 
			     (double *) w);
      break;
    case GS_OP_MUL:
      cudaMemset(v, 1, (*m) * sizeof(double));
      gather_kernel_mul<double>
	<<<nblcks, nthrds>>>((double *) v, *m, *o, (int *) dg,
			     (double *) u, *n, (int *) gd, 
			     (double *) w);
      break;
    case GS_OP_MIN:
      cudaMemset(v, INT_MAX, (*m) * sizeof(double));
      gather_kernel_min<double>
	<<<nblcks, nthrds>>>((double *) v, *m, *o, (int *) dg,
			     (double *) u, *n, (int *) gd, 
			     (double *) w);
      break;
    case GS_OP_MAX:
      cudaMemset(v, -INT_MAX, (*m) * sizeof(double));
      gather_kernel_max<double>
	<<<nblcks, nthrds>>>((double *) v, *m, *o, (int *) dg,
			     (double *) u, *n, (int *) gd, 
			     (double *) w);
      break;
    }
  }

  /**
   * Fortran wrapper for device scatter kernel
   */
  void cuda_scatter_kernel(void *v, int *m, void *dg,
			  void *u, int *n, void *gd, void *w) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    scatter_kernel<double>
      <<<nblcks, nthrds>>>((double *) v, *m, (int *) dg,
			   (double *) u, *n, (int *) gd, 
			   (double *) w);
  }
}
