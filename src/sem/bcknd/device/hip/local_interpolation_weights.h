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
#ifndef __HIP_LOCAL_INTERPOLATION_WEIGHTS_H__
#define __HIP_LOCAL_INTERPOLATION_WEIGHTS_H__

#include <device/device_config.h>

// Fornberg's recurrence for derivative order zero, matching fd_weights_full.
template<int LX>
__device__ inline void local_interpolation_weights_1d(
    real xi, const real *nodes, real *weights, int point) {
  real_xp c[LX];
  for (int j = 0; j < LX; ++j) c[j] = real_xp(0);
  c[0] = real_xp(1);
  real_xp c1 = real_xp(1);
  real_xp c4 = real_xp(nodes[0] - xi);
  for (int i = 1; i < LX; ++i) {
    real_xp c2 = real_xp(1);
    const real_xp c5 = c4;
    c4 = real_xp(nodes[i] - xi);
    for (int j = 0; j < i; ++j) {
      const real_xp c3 = real_xp(nodes[i] - nodes[j]);
      c2 *= c3;
      c[i] = -c1*c5*c[i-1]/c2;
      c[j] = c4*c[j]/c3;
    }
    c1 = c2;
  }
  for (int j = 0; j < LX; ++j) weights[point*LX+j] = real(c[j]);
}

template<int LX>
__device__ inline void local_interpolation_compute_weights_point(
    real r, real s, real t, const real *zg, real *wr, real *ws,
    real *wt, int p) {
  if (r <= real(1.1) && r >= real(-1.1) &&
      s <= real(1.1) && s >= real(-1.1) &&
      t <= real(1.1) && t >= real(-1.1)) {
    local_interpolation_weights_1d<LX>(r, zg, wr, p);
    local_interpolation_weights_1d<LX>(s, zg+LX, ws, p);
    local_interpolation_weights_1d<LX>(t, zg+2*LX, wt, p);
  } else {
    for (int j = 0; j < LX; ++j) {
      wr[p*LX+j] = real(0);
      ws[p*LX+j] = real(0);
      wt[p*LX+j] = real(0);
    }
  }
}

template<int LX>
__global__ void local_interpolation_compute_weights_kernel(
    const real *rst, const real *zg, real *wr, real *ws, real *wt, int n) {
  const int p = blockIdx.x*blockDim.x + threadIdx.x;
  if (p >= n) return;
  local_interpolation_compute_weights_point<LX>(
      rst[3*p], rst[3*p+1], rst[3*p+2], zg, wr, ws, wt, p);
}

template<int LX>
__global__ void local_interpolation_compute_weights_3arrays_kernel(
    const real *r, const real *s, const real *t, const real *zg,
    real *wr, real *ws, real *wt, int n) {
  const int p = blockIdx.x*blockDim.x + threadIdx.x;
  if (p >= n) return;
  local_interpolation_compute_weights_point<LX>(
      r[p], s[p], t[p], zg, wr, ws, wt, p);
}

#endif
