#include <hip/hip_runtime.h>
#include <cstdio>

#define GS_OP_ADD  1
#define GS_OP_MUL  2
#define GS_OP_MIN  3
#define GS_OP_MAX  4

/**
 * Device gather kernel for addition of data
 * \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
 */
__global__ void gather_kernel_add(double *v, int m, int o, int *dg,
				  double *u, int n, int *gd, double *w) {
}

/**
 * Device gather kernel for multiplication of data
 * \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
 */
__global__ void gather_kernel_mul(double *v, int m, int o, int *dg,
				  double *u, int n, int *gd, double *w) {
}

/**
 * Device gather kernel for minimum of data
 * \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
 */
__global__ void gather_kernel_min(double *v, int m, int o, int *dg,
				  double *u, int n, int *gd, double *w) {
}

/**
 * Device gather kernel for maximum of data
 * \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
 */
__global__ void gather_kernel_max(double *v, int m, int o, int *dg,
				  double *u, int n, int *gd, double *w) {
}

/**
 * Device scatter kernel
 * \f$ u(gd(i) = v(dg(i)) \f$
 */
__global__ void scatter_kernel(double *v, int m, int *dg,
			       double *u, int n, int *gd, double *w) {
}

extern "C" {

  /** 
   * Fortran wrapper for device gather kernels
   */
  void hip_gather_kernel(void *v, int *m, int *o, void *dg,
			 void *u, int *n, void *gd, void *w, int *op) {

    const int nblcks = (*m) / 256;
    const int nthrds = 256;

    switch (*op) {
    case GS_OP_ADD:
      hipLaunchKernelGGL(gather_kernel_add,
			 dim3(nblcks), dim3(nthrds), 0, 0,
			 (double *) v, *m, *o, (int *) dg,
			 (double *) u, *n, (int *) gd, (double *) w);
      break;
    case GS_OP_MUL:
      hipLaunchKernelGGL(gather_kernel_mul,
			 dim3(nblcks), dim3(nthrds), 0, 0,
			 (double *) v, *m, *o, (int *) dg,
			 (double *) u, *n, (int *) gd, (double *) w);
      break;
    case GS_OP_MIN:
      hipLaunchKernelGGL(gather_kernel_min,
			 dim3(nblcks), dim3(nthrds), 0, 0,
			 (double *) v, *m, *o, (int *) dg,
			 (double *) u, *n, (int *) gd, (double *) w);
      break;
    case GS_OP_MAX:
      hipLaunchKernelGGL(gather_kernel_max,
			 dim3(nblcks), dim3(nthrds), 0, 0,
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

    const int nblcks = (*m) / 256;
    const int nthrds = 256;

    hipLaunchKernelGGL(scatter_kernel,
		       dim3(nblcks), dim3(nthrds), 0, 0,
		       (double *) v, *m, (int *) dg,
		       (double *) u, *n, (int *) gd, (double *) w);
  }
}
