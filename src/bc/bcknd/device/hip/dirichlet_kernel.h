
/**
 * Device kernel for scalar apply for a Dirichlet condition
 */
template< typename T >
__global__ void dirichlet_apply_scalar_kernel(const int * __restrict__ msk,
					      T * __restrict__ x,
					      const T g,
					      const int m) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = msk[i] -1;
    x[k] = g;
  }
}

/**
 * Device kernel for vector apply for a Dirichlet condition
 */
template< typename T >
__global__ void dirichlet_apply_vector_kernel(const int * __restrict__ msk,
					      T * __restrict__ x,
					      T * __restrict__ y,
					      T * __restrict__ z,
					      const T g,
					      const int m) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = msk[i] -1;
    x[k] = g;
    y[k] = g;
    z[k] = g;
  }
}

