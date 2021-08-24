#include <climits>
#include <cstdio>

#define GS_OP_ADD  1
#define GS_OP_MUL  2
#define GS_OP_MIN  3
#define GS_OP_MAX  4

/**
 * Device gather kernel for addition of data
 * \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
 */
__global__ void gather_kernel_add(double * __restrict__ v,
				  const int m,
				  const int o,
				  const int * __restrict__ dg,
				  const double * __restrict__ u,
				  const int n,
				  const int * __restrict__ gd,
				  double * __restrict__ w) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < (abs(o) - 1); i += str) {
    w[i] = u[gd[i] - 1];
  }

  for (int i = idx; i < (abs(o) - 1); i += str) {
    atomicAdd(&v[dg[i] - 1], w[i]);
  }

  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	double tmp = u[gd[i] - 1] + u[gd[i+1] - 1];
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device gather kernel for multiplication of data
 * \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
 */
__global__ void gather_kernel_mul(double * __restrict__ v,
				  const int m,
				  const int o,
				  const int * __restrict__ dg,
				  const double * __restrict__ u,
				  const int n,
				  const int * __restrict__ gd,
				  double * __restrict__ w) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < (abs(o) - 1); i += str) {
    w[i] = u[gd[i] - 1];
  }

  for (int i = idx; i < (abs(o) - 1); i += str) {
    /** @todo Fix atomic mul
       atomicMul(&v[dg[i] - 1], w[i]);
    */
    v[dg[i] - 1] *= w[i];
  }

  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	double tmp = u[gd[i] - 1] * u[gd[i+1] - 1];
	v[dg[i] - 1] = tmp;
      }
    }
  }

}

/**
 * Device gather kernel for minimum of data
 * \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
 */
__global__ void gather_kernel_min(double * __restrict__ v,
				  const int m,
				  const int o,
				  const int * __restrict__ dg,
				  const double * __restrict__ u,
				  const int n,
				  const int * __restrict__ gd,
				  double * __restrict__ w) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < (abs(o) - 1); i += str) {
    w[i] = u[gd[i] - 1];
  }

  for (int i = idx; i < (abs(o) - 1); i += str) {
    /** @todo Fix atomic min
       atomicMin(&v[dg[i] - 1], w[i]);
    */
    v[dg[i] - 1] = min(v[dg[i] - 1], w[i]);
  }

  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	double tmp = min(u[gd[i] - 1], u[gd[i+1] - 1]);
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device gather kernel for maximum of data
 * \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
 */
__global__ void gather_kernel_max(double * __restrict__ v,
				  const int m,
				  const int o,
				  const int * __restrict__ dg,
				  const double * __restrict__ u,
				  const int n,
				  const int * __restrict__ gd,
				  double * __restrict__ w) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < (abs(o) - 1); i += str) {
    w[i] = u[gd[i] - 1];
  }

  for (int i = idx; i < (abs(o) - 1); i += str) {
    /** @todo Fix atomic max
       atomicMax(&v[dg[i] - 1], w[i]);
    */
    v[dg[i] - 1] = max(v[dg[i] - 1], w[i]);
  }

  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	double tmp = max(u[gd[i] - 1], u[gd[i+1] - 1]);
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device scatter kernel
 * \f$ u(gd(i) = v(dg(i)) \f$
 */
__global__ void scatter_kernel(double * __restrict__ v,
			       const int m,
			       const int * __restrict__ dg,
			       double * __restrict__ u,
			       const int n,
			       const int * __restrict__ gd,
			       double * __restrict__ w) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < m; i += str) {
    w[i] = v[dg[i] - 1];
  }

  for (int i = idx; i < m; i += str) {
    u[gd[i] - 1] = w[i];
  }
  
}

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
      gather_kernel_add<<<nblcks, nthrds>>>((double *) v, *m, *o, (int *) dg,
					    (double *) u, *n, (int *) gd, 
					    (double *) w);
      break;
    case GS_OP_MUL:
      cudaMemset(v, 1, (*m) * sizeof(double));
      gather_kernel_mul<<<nblcks, nthrds>>>((double *) v, *m, *o, (int *) dg,
					    (double *) u, *n, (int *) gd, 
					    (double *) w);
      break;
    case GS_OP_MIN:
      cudaMemset(v, INT_MAX, (*m) * sizeof(double));
      gather_kernel_min<<<nblcks, nthrds>>>((double *) v, *m, *o, (int *) dg,
					    (double *) u, *n, (int *) gd, 
					    (double *) w);
      break;
    case GS_OP_MAX:
      cudaMemset(v, -INT_MAX, (*m) * sizeof(double));
      gather_kernel_max<<<nblcks, nthrds>>>((double *) v, *m, *o, (int *) dg,
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

    scatter_kernel<<<nblcks, nthrds>>>((double *) v, *m, (int *) dg,
				       (double *) u, *n, (int *) gd, 
				       (double *) w);
  }
}
