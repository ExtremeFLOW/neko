
/**
 * Device kernel for add2s1
 */
template< typename T >
__global__ void add2s1_kernel(T * __restrict__ a,
			      const T * __restrict__ b,
			      const T c1,
			      const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = c1 * a[i] + b[i];
  }
}

/**
 * Device kernel for add2s2
 */
template< typename T >
__global__ void add2s2_kernel(T * __restrict__ a,
			      const T * __restrict__ b,
			      const T c1,
			      const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c1 * b[i];
  }
}

/**
 * Device kernel for invcol2
 */
template< typename T >
__global__ void invcol2_kernel(T * __restrict__ a,
			       const T * __restrict__ b,
			       const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    a[i] = a[i] / b[i];
  }  
}

/** 
 * Device kernel for col2
 */
template< typename T >
__global__ void col2_kernel(T * __restrict__ a,
			    const T * __restrict__ b,
			    const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] * b[i];
  }  
}

/** 
 * Device kernel for col3
 */
template< typename T >
__global__ void col3_kernel(T * __restrict__ a,
			    const T * __restrict__ b,
			    const T * __restrict__ c,
			    const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] * c[i];
  }  
}

/** 
 * Device kernel for sub3
 */
template< typename T >
__global__ void sub3_kernel(T * __restrict__ a,
			    const T * __restrict__ b,
			    const T * __restrict__ c,
			    const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] - c[i];
  }  
}

/**
 * Device kernel for addcol3
 */
template< typename T >
__global__ void addcol3_kernel(T * __restrict__ a,
			    const T * __restrict__ b,
			    const T * __restrict__ c,
			    const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + b[i] * c[i];
  }  
  

}

/**
 * Device kernel for glsc3
 */
template< typename T >
__global__ void glsc3_kernel(const T * a,
			     const T * b,
			     const T * c,
			     T * buf_h,
			     const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  __shared__ T buf[1024];
  T tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    tmp += a[i] * b[i] * c[i];
  }
  buf[threadIdx.x] = tmp;
  __syncthreads();

  int i = blockDim.x>>1;
  while (i != 0) {
    if (threadIdx.x < i) {
      buf[threadIdx.x] += buf[threadIdx.x + i];
    }
    __syncthreads();
    i = i>>1;
  }
 
  if (threadIdx.x == 0) {
    buf_h[blockIdx.x] = buf[0];
  }
}

/**
 * Device kernel for glsc2
 */
template< typename T >
__global__ void glsc3_kernel(const T * a,
			     const T * b,
			     T * buf_h,
			     const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  __shared__ T buf[1024];
  T tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    tmp += a[i] * b[i];
  }
  buf[threadIdx.x] = tmp;
  __syncthreads();

  int i = blockDim.x>>1;
  while (i != 0) {
    if (threadIdx.x < i) {
      buf[threadIdx.x] += buf[threadIdx.x + i];
    }
    __syncthreads();
    i = i>>1;
  }
 
  if (threadIdx.x == 0) {
    buf_h[blockIdx.x] = buf[0];
  }
}

