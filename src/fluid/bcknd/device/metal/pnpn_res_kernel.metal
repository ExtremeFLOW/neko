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
 * Metal compute kernels for the incompressible Pn/Pn residuals.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Pressure residual, part 1.
 */
kernel void prs_res_part1_kernel(device float *ta1[[ buffer(0) ]],
                                 device float *ta2[[ buffer(1) ]],
                                 device float *ta3[[ buffer(2) ]],
                                 device const float *wa1[[ buffer(3) ]],
                                 device const float *wa2[[ buffer(4) ]],
                                 device const float *wa3[[ buffer(5) ]],
                                 device const float *f_u[[ buffer(6) ]],
                                 device const float *f_v[[ buffer(7) ]],
                                 device const float *f_w[[ buffer(8) ]],
                                 device const float *B[[ buffer(9) ]],
                                 device float *h1[[ buffer(10) ]],
                                 constant float &mu[[ buffer(11) ]],
                                 constant float &rho[[ buffer(12) ]],
                                 constant int &n[[ buffer(13) ]],
                                 uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float inv_rho = 1.0f / rho;
  h1[idx] = inv_rho;
  ta1[idx] = (f_u[idx] / rho) - ((wa1[idx] * (mu / rho)) * B[idx]);
  ta2[idx] = (f_v[idx] / rho) - ((wa2[idx] * (mu / rho)) * B[idx]);
  ta3[idx] = (f_w[idx] / rho) - ((wa3[idx] * (mu / rho)) * B[idx]);
}

/**
 * Pressure residual, part 2.
 */
kernel void prs_res_part2_kernel(device float *p_res[[ buffer(0) ]],
                                 device const float *wa1[[ buffer(1) ]],
                                 device const float *wa2[[ buffer(2) ]],
                                 device const float *wa3[[ buffer(3) ]],
                                 constant int &n[[ buffer(4) ]],
                                 uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  p_res[idx] = (-p_res[idx]) + (wa1[idx] + wa2[idx] + wa3[idx]);
}

/**
 * Pressure residual, part 3.
 */
kernel void prs_res_part3_kernel(device float *p_res[[ buffer(0) ]],
                                 device const float *ta1[[ buffer(1) ]],
                                 device const float *ta2[[ buffer(2) ]],
                                 device const float *ta3[[ buffer(3) ]],
                                 constant float &dtbd[[ buffer(4) ]],
                                 constant int &n[[ buffer(5) ]],
                                 uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  p_res[idx] = p_res[idx] - (dtbd * (ta1[idx] + ta2[idx] + ta3[idx]));
}

/**
 * Velocity residual update.
 */
kernel void vel_res_update_kernel(device float *u_res[[ buffer(0) ]],
                                  device float *v_res[[ buffer(1) ]],
                                  device float *w_res[[ buffer(2) ]],
                                  device const float *ta1[[ buffer(3) ]],
                                  device const float *ta2[[ buffer(4) ]],
                                  device const float *ta3[[ buffer(5) ]],
                                  device const float *f_u[[ buffer(6) ]],
                                  device const float *f_v[[ buffer(7) ]],
                                  device const float *f_w[[ buffer(8) ]],
                                  constant int &n[[ buffer(9) ]],
                                  uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  u_res[idx] = (-u_res[idx]) - ta1[idx] + f_u[idx];
  v_res[idx] = (-v_res[idx]) - ta2[idx] + f_v[idx];
  w_res[idx] = (-w_res[idx]) - ta3[idx] + f_w[idx];
}
