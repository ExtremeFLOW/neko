
/**
 * Device kernel for vector apply for a symmetry condition
 */
template< typename T >
__global__ void symmetry_apply_vector_kernel(const int * __restrict__ xmsk,
					     const int * __restrict__ ymsk,
					     const int * __restrict__ zmsk,
					     T * __restrict__ x,
					     T * __restrict__ y,
					     T * __restrict__ z,
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


