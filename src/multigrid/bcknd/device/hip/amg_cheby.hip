/*
 Copyright (c) 2025, The Neko Authors
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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>

template< typename T >
__global__ void amg_cheby_solve_part1(T * __restrict__ r,
                                      T * __restrict__ f,
                                      T * __restrict__ w,
                                      T * __restrict__ x,
                                      T * __restrict__ d,
                                      const T inv_thet,
                                      const bool zero_initial,
                                      const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const real rt = f[i] - (!zero_initial ? w[i] : 0.0);
    const real dt = rt * inv_thet;
    r[i] = rt;
    d[i] = dt;
    x[i] = x[i] + dt;
    
  }
}

template< typename T >
__global__ void amg_cheby_solve_part2(T * __restrict__ r,
                                      T * __restrict__ w,
                                      T * __restrict__ d,
                                      T * __restrict__ x,
                                      const T tmp1,
                                      const T tmp2,
                                      const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const real rt = r[i] - w[i];
    const real dt = tmp1 * d[i] + tmp2 * rt;
    r[i] = rt;
    d[i] = dt;
    x[i] = x[i] + dt;
  }
}

extern "C" {

  /**
   * Fortran wrapper for part1 of cheby solve
   */
  void hip_amg_cheby_solve_part1(void *r, void *f, void *w, void *x, void *d,
                                 real *inv_thet, int *n, bool *zero_initial,
                                 hipStream_t strm) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(amg_cheby_solve_part1<real>),
                       nblcks, nthrds, 0, strm, (real *) r, (real *) f,
                       (real *) w, (real *) x, (real *) d,
                       *inv_thet, *zero_initial, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for part2 of cheby solve
   */
  void hip_amg_cheby_solve_part2(void *r, void *w, void *d, void *x,
                                 real *tmp1, real *tmp2, int *n,
                                 hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(amg_cheby_solve_part2<real>),
                       nblcks, nthrds, 0, strm, (real *) r,
                       (real *) w, (real *) d, (real *) x, *tmp1, *tmp2, *n);
    HIP_CHECK(hipGetLastError());
  }

}
