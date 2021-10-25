
/**
 * Device gather kernel for addition of data
 * \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
 */
template< typename T >
__global__ void gather_kernel_add(T * __restrict__ v,
				  const int m,
				  const int o,
				  const int * __restrict__ dg,
				  const T * __restrict__ u,
				  const int n,
				  const int * __restrict__ gd,
				  T * __restrict__ w) {

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
	T tmp = u[gd[i] - 1] + u[gd[i+1] - 1];
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device gather kernel for multiplication of data
 * \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
 */
template< typename T >
__global__ void gather_kernel_mul(T * __restrict__ v,
				  const int m,
				  const int o,
				  const int * __restrict__ dg,
				  const T * __restrict__ u,
				  const int n,
				  const int * __restrict__ gd,
				  T * __restrict__ w) {

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
	T tmp = u[gd[i] - 1] * u[gd[i+1] - 1];
	v[dg[i] - 1] = tmp;
      }
    }
  }

}

/**
 * Device gather kernel for minimum of data
 * \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
 */
template< typename T >
__global__ void gather_kernel_min(T * __restrict__ v,
				  const int m,
				  const int o,
				  const int * __restrict__ dg,
				  const T * __restrict__ u,
				  const int n,
				  const int * __restrict__ gd,
				  T * __restrict__ w) {

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
	T tmp = min(u[gd[i] - 1], u[gd[i+1] - 1]);
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device gather kernel for maximum of data
 * \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
 */
template< typename T >
__global__ void gather_kernel_max(T * __restrict__ v,
				  const int m,
				  const int o,
				  const int * __restrict__ dg,
				  const T * __restrict__ u,
				  const int n,
				  const int * __restrict__ gd,
				  T * __restrict__ w) {

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
	T tmp = max(u[gd[i] - 1], u[gd[i+1] - 1]);
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device scatter kernel
 * \f$ u(gd(i) = v(dg(i)) \f$
 */
template< typename T >
__global__ void scatter_kernel(T * __restrict__ v,
			       const int m,
			       const int * __restrict__ dg,
			       T * __restrict__ u,
			       const int n,
			       const int * __restrict__ gd,
			       T * __restrict__ w) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < m; i += str) {
    w[i] = v[dg[i] - 1];
  }

  for (int i = idx; i < m; i += str) {
    u[gd[i] - 1] = w[i];
  }
  
}

