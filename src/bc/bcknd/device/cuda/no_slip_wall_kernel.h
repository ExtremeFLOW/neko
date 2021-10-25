/**
 * Device kernel for scalar apply for a no-slip wall conditon
 */
template< typename T >
__global__ void no_slip_wall_apply_scalar_kernel(const int * __restrict__ msk,
						 T * __restrict__ x,
						 const int m) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = 0.0;
  }
}

/**
 * Device kernel for vector apply for a no-slip wall conditon
 */
template< typename T >
__global__ void no_slip_wall_apply_vector_kernel(const int * __restrict__ msk,
						 T * __restrict__ x,
						 T * __restrict__ y,
						 T * __restrict__ z,
						 const int m) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = 0.0;
    y[k] = 0.0;
    z[k] = 0.0;
  }
}

