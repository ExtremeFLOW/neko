#ifndef __KRYLOV_PIPECG_KERNEL_H__
#define __KRYLOV_PIPECG_KERNEL_H__
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

#include <math/bcknd/device/hip/math_kernel.h>

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

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;
  
  __shared__ T buf1[64];
  __shared__ T buf2[64];
  __shared__ T buf3[64];
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

  tmp1 = reduce_warp<T>(tmp1);
  tmp2 = reduce_warp<T>(tmp2);
  tmp3 = reduce_warp<T>(tmp3);
  if (lane == 0) {
    buf1[wid] = tmp1;
    buf2[wid] = tmp2;
    buf3[wid] = tmp3;
  }
  __syncthreads();

  tmp1 = (threadIdx.x < blockDim.x / warpSize) ? buf1[lane] : 0;
  tmp2 = (threadIdx.x < blockDim.x / warpSize) ? buf2[lane] : 0;
  tmp3 = (threadIdx.x < blockDim.x / warpSize) ? buf3[lane] : 0;
  if (wid == 0) {
    tmp1 = reduce_warp<T>(tmp1);
    tmp2 = reduce_warp<T>(tmp2);
    tmp3 = reduce_warp<T>(tmp3);
  }

  if (threadIdx.x == 0) {
    buf_h1[blockIdx.x] = tmp1;
    buf_h2[blockIdx.x] = tmp2;
    buf_h3[blockIdx.x] = tmp3;
  }
  
}

#endif // __KRYLOV_PIPECG_KERNEL_H__
