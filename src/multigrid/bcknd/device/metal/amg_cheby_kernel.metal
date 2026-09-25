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

/**
 * Metal compute kernels for the tree-AMG Chebyshev coarse-grid solve.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Chebyshev solve, part 1.
 */
kernel void amg_cheby_solve_part1_kernel(device float *r[[ buffer(0) ]],
                                         device const float *f[[ buffer(1) ]],
                                         device const float *w[[ buffer(2) ]],
                                         device float *x[[ buffer(3) ]],
                                         device float *d[[ buffer(4) ]],
                                         constant float &inv_thet[[ buffer(5) ]],
                                         constant int &zero_initial[[ buffer(6) ]],
                                         constant int &n[[ buffer(7) ]],
                                         uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float rt = f[idx] - (zero_initial ? 0.0f : w[idx]);
  const float dt = rt * inv_thet;
  r[idx] = rt;
  d[idx] = dt;
  x[idx] = x[idx] + dt;
}

/**
 * Chebyshev solve, part 2.
 */
kernel void amg_cheby_solve_part2_kernel(device float *r[[ buffer(0) ]],
                                         device const float *w[[ buffer(1) ]],
                                         device float *d[[ buffer(2) ]],
                                         device float *x[[ buffer(3) ]],
                                         constant float &tmp1[[ buffer(4) ]],
                                         constant float &tmp2[[ buffer(5) ]],
                                         constant int &n[[ buffer(6) ]],
                                         uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float rt = r[idx] - w[idx];
  const float dt = tmp1 * d[idx] + tmp2 * rt;
  r[idx] = rt;
  d[idx] = dt;
  x[idx] = x[idx] + dt;
}
