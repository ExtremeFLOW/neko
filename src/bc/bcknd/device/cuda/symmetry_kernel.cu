#include <algorithm>

/**
 * Device kernel for vector apply for a symmetry condition
 */
__global__ void symmetry_apply_vector_kernel(const int * __restrict__ xmsk,
					     const int * __restrict__ ymsk,
					     const int * __restrict__ zmsk,
					     double * __restrict__ x,
					     double * __restrict__ y,
					     double * __restrict__ z,
					     const int m,
					     const int n,
					     const int l) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (xmsk[i] - 1);
    x[k] = 0.0;
  }

  for (int i = (idx + 1); i < n; i += str) {
    const int k = (ymsk[i] - 1);
    y[k] = 0.0;
  }

  for (int i = (idx + 1); i < l; i += str) {
    const int k = (zmsk[i] - 1);
    z[k] = 0.0;
  }
  
}

extern "C" {

  /** 
   * Fortran wrapper for device symmetry apply vector
   */
  void cuda_symmetry_apply_vector(void *xmsk, void *ymsk, void *zmsk,
				 void *x, void *y, void *z,
				 int *m, int *n, int *l) {

    const int max_len = std::max(std::max(*m, *n), *l);
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((max_len) + 1024 - 1)/ 1024, 1, 1);

    symmetry_apply_vector_kernel<<<nblcks, nthrds>>>((int *) xmsk,
						     (int *) ymsk,
						     (int *) zmsk,
						     (double *) x,
						     (double *) y,
						     (double *) z,
						     *m, *n, *l);
  }
 
}
