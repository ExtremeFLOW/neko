#include <hip/hip_runtime.h>
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
  void hip_gather_kernel(void *v, int *m, int *o, void *dg,
			 void *u, int *n, void *gd, void *w, int *op) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+ 1024 - 1)/ 1024, 1, 1);

    switch (*op) {
    case GS_OP_ADD:
      hipMemset(v, 0, (*m) * sizeof(double));
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_kernel_add<double>),
			 nblcks, nthrds, 0, 0,
			 (double *) v, *m, *o, (int *) dg,
			 (double *) u, *n, (int *) gd, (double *) w);
      break;
    case GS_OP_MUL:
      hipMemset(v, 1, (*m) * sizeof(double));
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_kernel_mul<double>),
			 nblcks, nthrds, 0, 0,
			 (double *) v, *m, *o, (int *) dg,
			 (double *) u, *n, (int *) gd, (double *) w);
      break;
    case GS_OP_MIN:
      hipMemset(v, INT_MAX, (*m) * sizeof(double));
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_kernel_min<double>),
			 nblcks, nthrds, 0, 0,
			 (double *) v, *m, *o, (int *) dg,
			 (double *) u, *n, (int *) gd, (double *) w);
      break;
    case GS_OP_MAX:
      hipMemset(v, -INT_MAX, (*m) * sizeof(double));
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_kernel_max<double>),
			 nblcks, nthrds, 0, 0,
			 (double *) v, *m, *o, (int *) dg,
			 (double *) u, *n, (int *) gd, (double *) w);
      break;
    }
  }

  /**
   * Fortran wrapper for device scatter kernel
   */
  void hip_scatter_kernel(void *v, int *m, void *dg,
			  void *u, int *n, void *gd, void *w) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(scatter_kernel<double>),
		       nblcks, nthrds, 0, 0,
		       (double *) v, *m, (int *) dg,
		       (double *) u, *n, (int *) gd, (double *) w);
  }
}
