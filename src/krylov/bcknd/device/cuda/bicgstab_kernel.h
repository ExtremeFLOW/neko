/*
 Copyright (c) 2026, The Neko Authors
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

#ifndef __KRYLOV_BICGSTAB_KERNEL_H__
#define __KRYLOV_BICGSTAB_KERNEL_H__

#include <math/bcknd/device/cuda/math_kernel.h>

/**
 * Kernel for the BiCGStab search direction update
 * \f$ p = r + \beta (p - \omega v) \f$
 */
template< typename T >
__global__ void bicgstab_update_p_kernel(T * __restrict__ p,
                                         const T * __restrict__ r,
                                         const T * __restrict__ v,
                                         const T beta,
                                         const T omega,
                                         const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    p[i] = r[i] + beta * (p[i] - omega * v[i]);
  }
}

/**
 * Kernel for a weighted inner product and squared norm in one pass,
 * \f$ (a^T M b, b^T M b) \f$
 *
 * The two partial sums of each block are stored interleaved, so that
 * glsc3_reduce_kernel with j = 2 leaves the results in buf_h[0..1].
 */
template< typename T, typename T_acc >
__global__ void bicgstab_product_and_norm_kernel(const T * __restrict__ a,
                                                 const T * __restrict__ b,
                                                 const T * __restrict__ mult,
                                                 T_acc * __restrict__ buf_h,
                                                 const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T_acc shared_prod[32];
  __shared__ T_acc shared_nrm[32];
  T_acc prod = 0.0;
  T_acc nrm = 0.0;
  for (int i = idx; i < n; i += str) {
    const T_acc w = static_cast<T_acc>(mult[i] * b[i]);
    prod += static_cast<T_acc>(a[i]) * w;
    nrm += static_cast<T_acc>(b[i]) * w;
  }

  prod = reduce_warp<T_acc>(prod);
  nrm = reduce_warp<T_acc>(nrm);
  if (lane == 0) {
    shared_prod[wid] = prod;
    shared_nrm[wid] = nrm;
  }
  __syncthreads();

  prod = (threadIdx.x < blockDim.x / warpSize) ? shared_prod[lane] : 0.0;
  nrm = (threadIdx.x < blockDim.x / warpSize) ? shared_nrm[lane] : 0.0;
  if (wid == 0) {
    prod = reduce_warp<T_acc>(prod);
    nrm = reduce_warp<T_acc>(nrm);
  }

  if (threadIdx.x == 0) {
    buf_h[blockIdx.x * 2] = prod;
    buf_h[blockIdx.x * 2 + 1] = nrm;
  }
}

/**
 * Kernel for BiCGStab part 1, \f$ s = r - \alpha v \f$,
 * reducing \f$ s^T M s \f$
 */
template< typename T, typename T_acc >
__global__ void bicgstab_part1_kernel(T * __restrict__ s,
                                      const T * __restrict__ r,
                                      const T * __restrict__ v,
                                      const T * __restrict__ mult,
                                      T_acc * __restrict__ buf_h,
                                      const T alpha,
                                      const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T_acc shared[32];
  T_acc sum = 0.0;
  for (int i = idx; i < n; i += str) {
    const T si = r[i] - alpha * v[i];
    s[i] = si;
    sum += static_cast<T_acc>(si * mult[i] * si);
  }

  sum = reduce_warp<T_acc>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0.0;
  if (wid == 0)
    sum = reduce_warp<T_acc>(sum);

  if (threadIdx.x == 0)
    buf_h[blockIdx.x] = sum;
}

/**
 * Kernel for BiCGStab part 2,
 * \f$ x = x + \alpha \hat{p} + \omega \hat{s} \f$ and
 * \f$ r = s - \omega t \f$, reducing \f$ (r^T M r, f^T M r) \f$
 *
 * The second value is the rho inner product of the next iteration, which
 * the recurrence would otherwise reduce separately at the top of the loop.
 */
template< typename T, typename T_acc >
__global__ void bicgstab_part2_kernel(T * __restrict__ x,
                                      T * __restrict__ r,
                                      const T * __restrict__ p_hat,
                                      const T * __restrict__ s_hat,
                                      const T * __restrict__ s,
                                      const T * __restrict__ t,
                                      const T * __restrict__ f,
                                      const T * __restrict__ mult,
                                      T_acc * __restrict__ buf_h,
                                      const T alpha,
                                      const T omega,
                                      const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T_acc shared_rtr[32];
  __shared__ T_acc shared_ftr[32];
  T_acc rtr = 0.0;
  T_acc ftr = 0.0;
  for (int i = idx; i < n; i += str) {
    x[i] = x[i] + alpha * p_hat[i] + omega * s_hat[i];
    const T ri = s[i] - omega * t[i];
    r[i] = ri;
    const T_acc w = static_cast<T_acc>(mult[i] * ri);
    rtr += static_cast<T_acc>(ri) * w;
    ftr += static_cast<T_acc>(f[i]) * w;
  }

  rtr = reduce_warp<T_acc>(rtr);
  ftr = reduce_warp<T_acc>(ftr);
  if (lane == 0) {
    shared_rtr[wid] = rtr;
    shared_ftr[wid] = ftr;
  }
  __syncthreads();

  rtr = (threadIdx.x < blockDim.x / warpSize) ? shared_rtr[lane] : 0.0;
  ftr = (threadIdx.x < blockDim.x / warpSize) ? shared_ftr[lane] : 0.0;
  if (wid == 0) {
    rtr = reduce_warp<T_acc>(rtr);
    ftr = reduce_warp<T_acc>(ftr);
  }

  if (threadIdx.x == 0) {
    buf_h[blockIdx.x * 2] = rtr;
    buf_h[blockIdx.x * 2 + 1] = ftr;
  }
}

#endif // __KRYLOV_BICGSTAB_KERNEL_H__
