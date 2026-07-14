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
 * Metal compute kernel for the compressible flow maximum wave speed.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Device kernel for computing the maximum wave speed |u| + c
 */
kernel void compute_max_wave_speed_kernel(device float *max_wave_speed[[ buffer(0) ]],
                                          device const float *u[[ buffer(1) ]],
                                          device const float *v[[ buffer(2) ]],
                                          device const float *w[[ buffer(3) ]],
                                          constant float &gamma[[ buffer(4) ]],
                                          device const float *p[[ buffer(5) ]],
                                          device const float *rho[[ buffer(6) ]],
                                          constant int &n[[ buffer(7) ]],
                                          uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float vel_mag = sqrt(u[idx] * u[idx]
                             + v[idx] * v[idx]
                             + w[idx] * w[idx]);
  const float sound_speed = sqrt(gamma * p[idx] / rho[idx]);
  max_wave_speed[idx] = vel_mag + sound_speed;
}
