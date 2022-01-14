/**
 * Kernel for back-substitution of x and update of p
 */
template< typename T >
__global__ void cg_update_xp_kernel(T  * __restrict__  x,
                                    T * __restrict__ p,
                                    T ** __restrict__ u,
                                    const T * alpha,
                                    const T * beta,
                                    const int p_cur,
                                    const int p_space,
                                    const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;


  for (int i = idx; i < n; i+= str) {
    T tmp = 0.0;
    int p_prev = p_space;
    for (int j = 0; j < p_cur; j ++) {
      p[i] = beta[j]*p[i] + u[p_prev][i];
      tmp += alpha[j]*p[i];
      p_prev = j;
    }
    x[i] += tmp;
    u[p_space][i] = u[p_space-1][i];
  }
}

/**
 * Device kernel for pipecg_vecops
 */
template< typename T >
__global__ void pipecg_vecops_kernel(T  * __restrict__  p,
                                     T * __restrict__ q,
                                     T * __restrict__ r,
                                     T * __restrict__ s,
                                     T * __restrict__ u1,
                                     T * __restrict__ u2,
                                     T * __restrict__ w,
                                     T * __restrict__ z,
                                     T * __restrict__ ni,
                                     T * __restrict__ mi,
                                     const T alpha,
                                     const T beta,
                                     const T * mult,
                                     T * buf_h1,
                                     T * buf_h2,
                                     T * buf_h3,
                                     const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  __shared__ T buf1[1024];
  __shared__ T buf2[1024];
  __shared__ T buf3[1024];
  T tmp1 = 0.0;
  T tmp2 = 0.0;
  T tmp3 = 0.0;

  for (int i = idx; i < n; i+= str) {
    z[i] = beta * z[i] + ni[i];
    q[i] = beta * q[i] + mi[i];
    s[i] = beta * s[i] + w[i];
    r[i] =  r[i] - alpha * s[i];
    u2[i] =  u1[i] - alpha * q[i];
    w[i] =  w[i] - alpha * z[i];
    tmp1 = tmp1 + r[i] * mult[i] * u2[i];
    tmp2 = tmp2 + w[i] * mult[i] * u2[i];
    tmp3 = tmp3 + r[i] * mult[i] * r[i];

  }
  buf1[threadIdx.x] = tmp1;
  buf2[threadIdx.x] = tmp2;
  buf3[threadIdx.x] = tmp3;
  __syncthreads();

  int i = blockDim.x>>1;
  while (i != 0) {
    if (threadIdx.x < i) {
      buf1[threadIdx.x] += buf1[threadIdx.x + i];
      buf2[threadIdx.x] += buf2[threadIdx.x + i];
      buf3[threadIdx.x] += buf3[threadIdx.x + i];
    }
    __syncthreads();
    i = i>>1;
  }
 
  if (threadIdx.x == 0) {
    buf_h1[blockIdx.x] = buf1[0];
    buf_h2[blockIdx.x] = buf2[0];
    buf_h3[blockIdx.x] = buf3[0];
  }
}
