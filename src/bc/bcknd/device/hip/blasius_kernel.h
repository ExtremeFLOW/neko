/**
 * Device kernel for vector apply for a Blasius profile
 */
template< typename T >
__global__ void blasius_apply_vector_kernel(const int * __restrict__ msk,
					    T * __restrict__ x,
					    T * __restrict__ y,
					    T * __restrict__ z,
					    const T * __restrict__ bla_x,
					    const T * __restrict__ bla_y,
					    const T * __restrict__ bla_z,
					    const int m) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < m; i += str) {
    const int k = msk[i + 1] - 1;
    x[k] = bla_x[i];
    y[k] = bla_y[i];
    z[k] = bla_z[i];
  }
}
