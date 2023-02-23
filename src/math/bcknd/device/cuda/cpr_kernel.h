/*
 Copyright (c) 2021-2022, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * Device kernel for lcsc3
 */
template< typename T, const int NXYZ >
__global__ void lcsc3_kernel(T * buf_h,
		             const T * a,
                             const T * b,
                             const T * c) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ T buf[NXYZ];
  T tmp = 0.0;

  tmp += a[idx] * b[idx] * c[idx];
  
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
 * Device kernel for lcsum
 */
template< typename T, const int NXYZ >
__global__ void lcsum_kernel(T * buf_h,
		           const T * a) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ T buf[NXYZ];
  T tmp = 0.0;

  tmp += a[idx];
  
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
 * Device kernel for lcmin
 */
template< typename T, const int NXYZ >
__global__ void lcmin_kernel(T * buf_h,
		             const T * a) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ T buf[NXYZ];
  T tmp = 0.0;
  T tmp1 = 0.0;
  T tmp2 = 0.0;

  tmp += a[idx];
  
  buf[threadIdx.x] = tmp;
  __syncthreads();

  int i = blockDim.x>>1;
  while (i != 0) {
    if (threadIdx.x < i) {
      tmp1 = buf[threadIdx.x];
      tmp2 = buf[threadIdx.x+i];
      if (tmp1>tmp2) {
        buf[threadIdx.x] = buf[threadIdx.x + i];
      }
    }
    __syncthreads();
    i = i>>1;
  }
 
  if (threadIdx.x == 0) {
    buf_h[blockIdx.x] = buf[0];
  }
}


/**
 * Device kernel for lcmax
 */
template< typename T, const int NXYZ >
__global__ void lcmax_kernel(T * buf_h,
		             const T * a) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ T buf[NXYZ];
  T tmp = 0.0;
  T tmp1 = 0.0;
  T tmp2 = 0.0;

  tmp += a[idx];
  
  buf[threadIdx.x] = tmp;
  __syncthreads();

  int i = blockDim.x>>1;
  while (i != 0) {
    if (threadIdx.x < i) {
      tmp1 = buf[threadIdx.x];
      tmp2 = buf[threadIdx.x+i];
      if (tmp1<tmp2) {
        buf[threadIdx.x] = buf[threadIdx.x + i];
      }
    }
    __syncthreads();
    i = i>>1;
  }
 
  if (threadIdx.x == 0) {
    buf_h[blockIdx.x] = buf[0];
  }
}


/**
 * Device kernel for lcsort_abs
 * Sort with absolute values
 * Using slow bubble sort with one thread sorting per block
 * buf_h and buf_k must be same size as a and key
 */
template< typename T, const int NXYZ >
__global__ void lcsort_abs_kernel(T * buf_h,
		                  int * buf_i,
		             const T * a,
		             const int * key) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ T buf[NXYZ];
  __shared__ int buf_k[NXYZ];
  T tmp = 0.0;
  T tmp1 = 0.0;
  T tmp2 = 0.0;
  T tmp3 = 0.0;
 
  buf[threadIdx.x]   = a[idx];
  buf_k[threadIdx.x] = key[idx];
  __syncthreads();

  if (threadIdx.x == 0) {
    for (int i = 0; i< NXYZ; i += 1) {

      for (int j = 0; j < (NXYZ-i-1); j += 1) {

      	tmp  = buf[j];
        tmp1 = buf[j+1];
        tmp2 = buf_k[j];
        tmp3 = buf_k[j+1];
        if (abs(tmp)<abs(tmp1)) {
          buf[j]     = tmp1;
          buf[j+1]   = tmp;
          buf_k[j]   = tmp3;
          buf_k[j+1] =tmp2;
        }
      } 
    }
  }  
  __syncthreads();

  buf_h[idx] = buf[threadIdx.x];
  buf_i[idx] = buf_k[threadIdx.x]; 
 
}


/**
 * Device kernel for lcsort_bykey
 * Put entries in the order given by the key contents
 * buf_h and buf_k must be same size as a and key
 */
template< typename T, const int NXYZ >
__global__ void lcsort_bykey_kernel(T * buf_h,
		                    int * buf_i,
				    const T * a,
		                    const int * key) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ T buf[NXYZ];
  __shared__ int buf_k[NXYZ];
  T tmp = 0.0;
  T tmp1 = 0.0;
  T tmp2 = 0.0;
  T tmp3 = 0.0;
 
  buf[threadIdx.x]   = a[idx];
  buf_k[threadIdx.x] = key[idx];
  __syncthreads();

  for (int i = 0; i< NXYZ; i += 1) {

      	tmp  = buf_k[i];

        if (threadIdx.x == tmp) {
          buf_h[idx]= buf[threadIdx.x];
          buf_i[idx]= buf_k[threadIdx.x];
        }
  }
  __syncthreads();

}

template< typename T >
__global__ void lctnsr3d_kernel(T  * __restrict__  v,
                              const int nv,
                              const T * __restrict__ u,
                              const int nu,
                              const T * __restrict__ A,
                              const T * __restrict__ Bt,
                              const T * __restrict__ Ct,
                              const T * __restrict__ bp,
                              const T * __restrict__ bp_key) {
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
      if (bp_key[e]>0.5){	    
        tmp += A[i+l*nv]*u[l+nu*j+e*nu*nu*nu];
      }
      else{
        tmp += bp[i+l*nv]*u[l+nu*j+e*nu*nu*nu];
      }
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
      if (bp_key[e]>0.5){	    
        tmp += Bt[l+j*nu]*shwork[l*nv+ik2];
      }
      else{
        tmp += bp[l+j*nu]*shwork[l*nv+ik2];
      }
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
      if (bp_key[e]>0.5){	    
        tmp += Ct[l+k*nu]*shwork2[ij2 + l*nv*nv];
      }
      else{
        tmp += bp[l+k*nu]*shwork2[ij2 + l*nv*nv];
      }
    }
    v[ijk+e*nv*nv*nv] = tmp;
  }
}





