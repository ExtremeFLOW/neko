template< typename T >
__global__ void tnsr3d_kernel(T  * __restrict__  v,
                                   const int nv,
                                   const T * __restrict__ u,
                                   const int nu,
                                   const T * __restrict__ A,
                                   const T * __restrict__ Bt,
                                   const T * __restrict__ Ct) {
  __shared__ T shwork[2048];
  __shared__ T shwork2[2048];
  const int idx = threadIdx.x;
  const int str = blockDim.x;
  const int e = blockIdx.x;
  for (int ii = idx; ii< nu*nu*nv; ii += str){
    T tmp = 0.0;
    int j = ii/nv;
    int i = ii - j*nv;
    for( int l = 0; l < nu; l++){
      tmp += A[i+l*nv]*u[l+nu*j+e*nu*nu*nu];
    }
    shwork[ii] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< nu*nv*nv; ijk += str){
    const int jk = ijk / nv;
    const int i = ijk - jk * nv;
    const int k = jk / nv;
    const int j = jk - k * nv;
    T tmp = 0.0;
    const int ik2 = i + k*nv*nu; 
    for( int l = 0; l < nu; l++){
      tmp += Bt[l+j*nu]*shwork[l*nv+ik2];
    }
    shwork2[ijk] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< nv*nv*nv; ijk += str){
    const int jk = ijk / nv;
    const int i = ijk - jk * nv;
    const int k = jk / nv;
    const int j = jk - k * nv;
    T tmp = 0.0;
    const int ij2 = i + j*nv; 
    for( int l = 0; l < nu; l++){
      tmp += Ct[l+k*nu]*shwork2[ij2 + l*nv*nv];
    }
    v[ijk+e*nv*nv*nv] = tmp;
  }
}


