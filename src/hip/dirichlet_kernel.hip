#include <hip/hip_runtime.h>

/**
 * Device kernel for scalar apply for a Dirichlet condition
 */
__global__ void dirichlet_apply_scalar_kernel(const int * __restrict__ msk,
					      double * __restrict__ x,
					      const double g,
					      const int m) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = g;
  }
}

/**
 * Device kernel for vector apply for a Dirichlet condition
 */
__global__ void dirichlet_apply_vector_kernel(const int * __restrict__ msk,
					      double * __restrict__ x,
					      double * __restrict__ y,
					      double * __restrict__ z,
					      const double g,
					      const int m) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = g;
    y[k] = g;
    z[k] = g;
  }
}

extern "C" {

  /** 
   * Fortran wrapper for device dirichlet apply scalar
   */
  void hip_dirichlet_apply_scalar(void *msk, void *x,
				  double *g, int *m) {
    
    const dim3 nthrds(512, 1, 1);
    const dim3 nblcks(((*m)+512 - 1)/ 512, 1, 1);

    hipLaunchKernelGGL(dirichlet_apply_scalar_kernel,
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (double *) x, *g, *m);
  }
  
  /** 
   * Fortran wrapper for device dirichlet apply vector
   */
  void hip_dirichlet_apply_vector(void *msk, void *x, void *y,
				  void *z, double *g, int *m) {
    
    const dim3 nthrds(512, 1, 1);
    const dim3 nblcks(((*m)+512 - 1)/ 512, 1, 1);

    hipLaunchKernelGGL(dirichlet_apply_vector_kernel,
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (double *) x,
		       (double *) y, (double *) z, *g, *m);
  }
 
}
