#ifndef __KRYLOV_GMRES_KERNEL_H__
#define __KRYLOV_GMRES_KERNEL_H__
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

#include <math/bcknd/device/cuda/math_kernel.h>

/**
 * Kernel for back-substitution of x and update of p
 */
template< typename T >
__global__ void gmres_part2_kernel(T  * __restrict__  w,
                                   T * const * __restrict__ v,
                                   const T * __restrict__ mult,
                                   const T * __restrict__ h,
                                   T * __restrict__ buf_h1,
                                   const int j,
                                   const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;
  
  __shared__ T shared[32];
  T tmp1 = 0.0;

  for (int i = idx; i < n; i+= str) {
    T tmp = 0.0;
    for (int k = 0; k < j; k ++) {
      tmp += -h[k]*v[k][i];
    }
    w[i] += tmp;
    tmp1 += w[i]*w[i]*mult[i];
  }

  tmp1 = reduce_warp<T>(tmp1);
  if (lane == 0)
    shared[wid] = tmp1;
  __syncthreads();

  tmp1 = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    tmp1 = reduce_warp<T>(tmp1);

  if (threadIdx.x == 0)
    buf_h1[blockIdx.x] = tmp1;

}



#endif // __KRYLOV_GMRES_KERNEL_H__
