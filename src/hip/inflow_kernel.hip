#include <hip/hip_runtime.h>

/**
 * Device kernel for vector apply for a Dirichlet condition
 */
__global__ void inflow_apply_vector_kernel(const int * __restrict__ msk,
					   double * __restrict__ x,
					   double * __restrict__ y,
					   double * __restrict__ z,
					   const double gx,
					   const double gy,
					   const double gz,
					   const int m) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = gx;
    y[k] = gy;
    z[k] = gz;
  }
}

extern "C" {

  /** 
   * Fortran wrapper for device inflow apply vector
   */
  void hip_inflow_apply_vector(void *msk, void *x, void *y,
				  void *z, void *g, int *m) {

    const double gx = ((double *)g)[0];
    const double gy = ((double *)g)[1];
    const double gz = ((double *)g)[2];
    
    const dim3 nthrds(512, 1, 1);
    const dim3 nblcks(((*m)+512 - 1)/ 512, 1, 1);

    hipLaunchKernelGGL(inflow_apply_vector_kernel,
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (double *) x,
		       (double *) y, (double *) z,
		       gx, gy, gz, *m);
  }
 
}
