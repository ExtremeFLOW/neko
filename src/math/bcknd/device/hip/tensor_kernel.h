#ifndef __MATH_TENSOR_KERNEL_H__
#define __MATH_TENSOR_KERNEL_H__
/*
 Copyright (c) 2022, The Neko Authors
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
template< typename T, const int N >
__global__ void tnsr3d_el_kernel(T  * __restrict__  v,
                              const int nv,
                              const T * __restrict__ u,
                              const int nu,
                              const T * __restrict__ A,
                              const T * __restrict__ Bt,
                              const T * __restrict__ Ct,
                              const int *  elements,
                              const int n_points) {
  __shared__ T shwork[N*N*N];
  __shared__ T shwork2[N*N*N];

  const int idx = threadIdx.x;
  const int str = blockDim.x;
  const int pt = blockIdx.x;
  const int e = elements[pt];

  for (int ii = idx; ii< nu*nu*nv; ii += str) {
    T tmp = 0.0;
    int j = ii/nv;
    int i = ii - j*nv;
    for( int l = 0; l < nu; l++){
      tmp += A[i+l*nv+pt*nv*nu]*u[l+nu*j+e*nu*nu*nu];
    }
    shwork[ii] = tmp;
  }

  __syncthreads();

  for (int ijk = idx; ijk< nu*nv*nv; ijk += str) {
    const int jk = ijk / nv;
    const int i = ijk - jk * nv;
    const int k = jk / nv;
    const int j = jk - k * nv;
    T tmp = 0.0;
    const int ik2 = i + k*nv*nu;
    for( int l = 0; l < nu; l++){
      tmp += Bt[l+j*nu+pt*nv*nu]*shwork[l*nv+ik2];
    }
    shwork2[ijk] = tmp;
  }

  __syncthreads();

  for (int ijk = idx; ijk< nv*nv*nv; ijk += str) {
    const int jk = ijk / nv;
    const int i = ijk - jk * nv;
    const int k = jk / nv;
    const int j = jk - k * nv;
    T tmp = 0.0;
    const int ij2 = i + j*nv;
    for( int l = 0; l < nu; l++){
      tmp += Ct[l+k*nu+pt*nv*nu]*shwork2[ij2 + l*nv*nv];
    }
    v[ijk+pt*nv*nv*nv] = tmp;
  }

}


template< typename T, const int N >
__global__ void tnsr3d_kernel(T  * __restrict__  v,
	                      const int nv,
                              const T * __restrict__ u,
	                      const int nu,
                              const T * __restrict__ A,
                              const T * __restrict__ Bt,
                              const T * __restrict__ Ct) {
  __shared__ T shwork[N*N*N];
  __shared__ T shwork2[N*N*N];
  
  const int idx = threadIdx.x;
  const int str = blockDim.x;
  const int e = blockIdx.x;
  
  for (int ii = idx; ii< nu*nu*nv; ii += str) {
    T tmp = 0.0;
    int j = ii/nv;
    int i = ii - j*nv;
    for( int l = 0; l < nu; l++){
      tmp += A[i+l*nv]*u[l+nu*j+e*nu*nu*nu];
    }
    shwork[ii] = tmp;
  }
  
  __syncthreads();
  
  for (int ijk = idx; ijk< nu*nv*nv; ijk += str) {
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
  
  for (int ijk = idx; ijk< nv*nv*nv; ijk += str) {
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



#endif // __MATH_TENSOR_KERNEL_H__
