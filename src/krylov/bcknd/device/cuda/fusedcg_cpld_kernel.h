#ifndef __KRYLOV_FUSEDCG_CPLD_KERNEL_H__
#define __KRYLOV_FUSEDCG_CPLD_KERNEL_H__
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

#include <math/bcknd/device/cuda/math_kernel.h>

/**
 * Kernel for first coupled fusedcg part
 */
template< typename T >
__global__ void fusedcg_cpld_part1_kernel(T * __restrict__  a1,
                                          T * __restrict__  a2,
                                          T * __restrict__  a3,
                                          T * __restrict__  b1,
                                          T * __restrict__  b2,
                                          T * __restrict__  b3,
                                          T * __restrict__  tmp,
                                          const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i+= str) {
    tmp[i] = a1[i]*b1[i] + a2[i]*b2[i] + a3[i]*b3[i];
  }

}

/**
 * Kernel for update of px
 */
template< typename T >
__global__ void fusedcg_cpld_update_p_kernel(T * __restrict__  p1,
                                             T * __restrict__  p2,
                                             T * __restrict__  p3,
                                             const T * __restrict__ z1,
                                             const T * __restrict__ z2,
                                             const T * __restrict__ z3,
                                             const T * __restrict__ po1,
                                             const T * __restrict__ po2,
                                             const T * __restrict__ po3,
                                             const T beta,
                                             const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i+= str) {
    p1[i] = beta*po1[i] + z1[i];
    p2[i] = beta*po2[i] + z2[i];
    p3[i] = beta*po3[i] + z3[i];
  }
  
}

/**
 * Kernel for update of xx
 */
template< typename T >
__global__ void fusedcg_cpld_update_x_kernel(T * __restrict__  x1,
                                             T * __restrict__  x2,
                                             T * __restrict__  x3,
                                             const T ** p1,
                                             const T ** p2,
                                             const T ** p3,
                                             const T * __restrict__ alpha,
                                             const int p_cur,
                                             const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i+= str) {
    T tmp1 = 0.0;
    T tmp2 = 0.0;
    T tmp3 = 0.0;
    for (int j = 0; j < p_cur; j ++) {
      tmp1 += p1[j][i] * alpha[j];
      tmp2 += p2[j][i] * alpha[j];
      tmp3 += p3[j][i] * alpha[j];
    }
    x1[i] += tmp1;
    x2[i] += tmp2;
    x3[i] += tmp3;
  }
  
}

/**
 * Device kernel for fusedcg_cpld_part2
 */
template< typename T>
__global__ void fusedcg_cpld_part2_kernel(T * __restrict__ a1,
                                          T * __restrict__ a2,
                                          T * __restrict__ a3,
                                          const T * __restrict__ b,
                                          const T * __restrict__ c1,
                                          const T * __restrict__ c2,
                                          const T * __restrict__ c3,
                                          const T alpha,
                                          T * buf_h,
                                          const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;
  
  __shared__ T buf[32];
  T tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    T rt1 = a1[i] - alpha * c1[i];
    T rt2 = a2[i] - alpha * c2[i];
    T rt3 = a3[i] - alpha * c3[i];
    tmp = tmp + ((rt1*rt1 + rt2*rt2 + rt3*rt3) * b[i]);
    a1[i] = rt1;
    a2[i] = rt2;
    a3[i] = rt3;
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

#endif // __KRYLOV_FUSEDCG_CPLD_KERNEL_H__
