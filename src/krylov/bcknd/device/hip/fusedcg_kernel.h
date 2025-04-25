#ifndef __KRYLOV_FUSEDCG_KERNEL_H__
#define __KRYLOV_FUSEDCG_KERNEL_H__
/*
 Copyright (c) 2021-2024, The Neko Authors
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
 * Kernel for update of p
 */
template< typename T >
__global__ void fusedcg_update_p_kernel(T * __restrict__  p,
                                    const T * __restrict__ z,
                                    const T * __restrict__ po,
                                    const T beta,
                                    const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i+= str) {
    p[i] = beta*po[i] + z[i];
  }
  
}

/**
 * Kernel for update of x
 */
template< typename T >
__global__ void fusedcg_update_x_kernel(T * __restrict__  x,
                                        const T ** p,
                                        const T * __restrict__ alpha,
                                        const int p_cur,
                                        const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i+= str) {
    T tmp = 0.0;
    for (int j = 0; j < p_cur; j ++) {
      tmp += p[j][i] * alpha[j];
    }
    x[i] += tmp;
  }
  
}

/**
 * Device kernel for fusedcg_part2
 */
template< typename T>
__global__ void fusedcg_part2_kernel(T * __restrict__ a,
                                     const T * __restrict__ b,
                                     const T * __restrict__ c,
                                     const T alpha,
                                     T * buf_h,
                                     const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;
  
  __shared__ T buf[64];
  T tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    T rt = a[i] - alpha * c[i];
    tmp = tmp + rt * b[i] * rt;
    a[i] = rt;
  }

  tmp = reduce_warp<T>(tmp);
  if (lane == 0) {
    buf[wid] = tmp;
  }
  __syncthreads();
  
  tmp = (threadIdx.x < blockDim.x / warpSize) ? buf[lane] : 0;
  if (wid == 0) {
    tmp = reduce_warp<T>(tmp);
  }

  if (threadIdx.x == 0) {
    buf_h[blockIdx.x] = tmp;
  }
}
                                     

#endif // __KRYLOV_FUSEDCG_KERNEL_H__
