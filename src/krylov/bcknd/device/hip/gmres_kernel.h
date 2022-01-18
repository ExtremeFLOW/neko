/**
 * Kernel for back-substitution of x and update of p
 */
template< typename T >
__global__ void gmres_part2_kernel(T  * __restrict__  w,
                                   T ** __restrict__ v,
                                   T * const mult,
                                   T * const h,
                                   T * buf_h1,
                                   const int j,
                                   const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  __shared__ T buf1[1024];
  T tmp1 = 0.0;

  for (int i = idx; i < n; i+= str) {
    T tmp = 0.0;
    for (int k = 0; k < j; k ++) {
      tmp += -h[k]*v[k][i];
    }
    w[i] += tmp;
    tmp1 += w[i]*w[i]*mult[i];
  }
  buf1[threadIdx.x] = tmp1;
  __syncthreads();

  int i = blockDim.x>>1;
  while (i != 0) {
    if (threadIdx.x < i) {
      buf1[threadIdx.x] += buf1[threadIdx.x + i];
    }
    __syncthreads();
    i = i>>1;
  }
 
  if (threadIdx.x == 0) {
    buf_h1[blockIdx.x] = buf1[0];
  }
}


