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
 * Metal compute kernels for the design-field mapping functions
 * (smooth step, step, inverse permeability) used by the Brinkman term.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Second order smooth step function.
 * t = clamp((x - edge0) / (edge1 - edge0), 0, 1);
 * x = t^3 * (t * (6 t - 15) + 10).
 */
kernel void smooth_step_kernel(device float *x[[ buffer(0) ]],
                               constant float &edge0[[ buffer(1) ]],
                               constant float &edge1[[ buffer(2) ]],
                               constant int &n[[ buffer(3) ]],
                               uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  float t = (x[idx] - edge0) / (edge1 - edge0);
  t = clamp(t, 0.0f, 1.0f);
  x[idx] = t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

/**
 * Step function: x = (x < edge) ? left : right.
 */
kernel void step_kernel(device float *x[[ buffer(0) ]],
                        constant float &edge[[ buffer(1) ]],
                        constant float &left[[ buffer(2) ]],
                        constant float &right[[ buffer(3) ]],
                        constant int &n[[ buffer(4) ]],
                        uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float x_i = x[idx];
  x[idx] = (x_i < edge) ? left : right;
}

/**
 * Inverse permeability mapping.
 * x = k_0 + (k_1 - k_0) * x * (q + 1) / (q + x).
 */
kernel void permeability_kernel(device float *x[[ buffer(0) ]],
                                constant float &k_0[[ buffer(1) ]],
                                constant float &k_1[[ buffer(2) ]],
                                constant float &q[[ buffer(3) ]],
                                constant int &n[[ buffer(4) ]],
                                uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float x_i = x[idx];
  x[idx] = k_0 + (k_1 - k_0) * x_i * (q + 1.0f) / (q + x_i);
}
