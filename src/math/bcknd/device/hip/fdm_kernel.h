template< typename T >
__global__ void fdm_do_fast_kernel(T  * __restrict__ e,
                                   T * __restrict__ r,
                                   T * __restrict__ s,
                                   T * __restrict__ d,
                                   const int nl) {
  __shared__ T shwork[2048];
  __shared__ T shwork2[2048];
  __shared__ T A[256];
  __shared__ T Bt[256];
  __shared__ T Ct[256];

  const int idx = threadIdx.x;
  const int str = blockDim.x;
  const int el = blockIdx.x;
  if( idx < nl*nl){
     A[idx] = s[idx+nl*nl+el*nl*nl*3*2];
    Bt[idx] = s[idx+2*nl*nl+el*nl*nl*3*2];
    Ct[idx] = s[idx+2*2*nl*nl+el*nl*nl*3*2];
  }    
  __syncthreads();

  for (int ii = idx; ii< nl*nl*nl; ii += str){
    T tmp = 0.0;
    int j = ii/nl;
    int i = ii - j*nl;
    for( int l = 0; l < nl; l++){
      tmp += A[i+l*nl]*r[l+nl*j+el*nl*nl*nl];
    }
    shwork[ii] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< nl*nl*nl; ijk += str){
    const int jk = ijk / nl;
    const int i = ijk - jk * nl;
    const int k = jk / nl;
    const int j = jk - k * nl;
    T tmp = 0.0;
    const int ik2 = i + k*nl*nl; 
    for( int l = 0; l < nl; l++){
      tmp += Bt[l+j*nl]*shwork[l*nl+ik2];
    }
    shwork2[ijk] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< nl*nl*nl; ijk += str){
    const int jk = ijk / nl;
    const int i = ijk - jk * nl;
    const int k = jk / nl;
    const int j = jk - k * nl;
    T tmp = 0.0;
    const int ij2 = i + j*nl; 
    for( int l = 0; l < nl; l++){
      tmp += Ct[l+k*nl]*shwork2[ij2 + l*nl*nl];
    }
    r[ijk+el*nl*nl*nl] = tmp*d[ijk+el*nl*nl*nl];
  }
  __syncthreads();
  if( idx < nl*nl){
     A[idx] = s[idx+el*nl*nl*3*2];
    Bt[idx] = s[idx+nl*nl+2*nl*nl+el*nl*nl*3*2];
    Ct[idx] = s[idx+nl*nl+2*2*nl*nl+el*nl*nl*3*2];
  }  
  __syncthreads();

  for (int ii = idx; ii< nl*nl*nl; ii += str){
    T tmp = 0.0;
    int j = ii/nl;
    int i = ii - j*nl;
    for( int l = 0; l < nl; l++){
      tmp += A[i+l*nl]*r[l+nl*j+el*nl*nl*nl];
    }
    shwork[ii] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< nl*nl*nl; ijk += str){
    const int jk = ijk / nl;
    const int i = ijk - jk * nl;
    const int k = jk / nl;
    const int j = jk - k * nl;
    T tmp = 0.0;
    const int ik2 = i + k*nl*nl; 
    for( int l = 0; l < nl; l++){
      tmp += Bt[l+j*nl]*shwork[l*nl+ik2];
    }
    shwork2[ijk] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< nl*nl*nl; ijk += str){
    const int jk = ijk / nl;
    const int i = ijk - jk * nl;
    const int k = jk / nl;
    const int j = jk - k * nl;
    T tmp = 0.0;
    const int ij2 = i + j*nl; 
    for( int l = 0; l < nl; l++){
      tmp += Ct[l+k*nl]*shwork2[ij2 + l*nl*nl];
    }
    e[ijk+el*nl*nl*nl] = tmp;
  }

}


