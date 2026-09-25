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
 * Metal compute kernels for entropy viscosity regularization.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

#define EV_MAX_THREADS 256

/**
 * Device kernel for computing the entropy residual (BDF time derivative)
 */
kernel void entropy_visc_compute_residual_kernel(device float *entropy_residual[[ buffer(0) ]],
                                                 device const float *S[[ buffer(1) ]],
                                                 device const float *S_lag1[[ buffer(2) ]],
                                                 device const float *S_lag2[[ buffer(3) ]],
                                                 device const float *S_lag3[[ buffer(4) ]],
                                                 constant float &bdf1[[ buffer(5) ]],
                                                 constant float &bdf2[[ buffer(6) ]],
                                                 constant float &bdf3[[ buffer(7) ]],
                                                 constant float &bdf4[[ buffer(8) ]],
                                                 constant float &dt[[ buffer(9) ]],
                                                 constant int &n[[ buffer(10) ]],
                                                 uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  entropy_residual[idx] = (bdf1 * S[idx]
                           - bdf2 * S_lag1[idx]
                           - bdf3 * S_lag2[idx]
                           - bdf4 * S_lag3[idx]) / dt;
}

/**
 * Device kernel for computing viscosity from the entropy residual
 */
kernel void entropy_visc_compute_viscosity_kernel(device float *reg_coeff[[ buffer(0) ]],
                                                  device const float *entropy_residual[[ buffer(1) ]],
                                                  device const float *h[[ buffer(2) ]],
                                                  constant float &c_avisc_entropy[[ buffer(3) ]],
                                                  constant float &n_S[[ buffer(4) ]],
                                                  constant int &n[[ buffer(5) ]],
                                                  uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  reg_coeff[idx] = c_avisc_entropy * h[idx] * h[idx]
                   * entropy_residual[idx] / n_S;
}

/**
 * Device kernel for applying the element-wise maximum,
 * one element per threadgroup
 */
kernel void entropy_visc_apply_element_max_kernel(device float *reg_coeff[[ buffer(0) ]],
                                                  constant int &lx3[[ buffer(1) ]],
                                                  uint el [[ threadgroup_position_in_grid ]],
                                                  uint tid [[ thread_index_in_threadgroup ]]) {
  threadgroup float sdata[EV_MAX_THREADS];

  const int offset = int(el) * lx3;

  float max_val = -1e30f;
  for (int i = int(tid); i < lx3; i += EV_MAX_THREADS) {
    max_val = fmax(max_val, reg_coeff[offset + i]);
  }

  sdata[tid] = max_val;
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int s = EV_MAX_THREADS / 2; s > 0; s >>= 1) {
    if (int(tid) < s) {
      sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }

  max_val = sdata[0];

  for (int i = int(tid); i < lx3; i += EV_MAX_THREADS) {
    reg_coeff[offset + i] = max_val;
  }
}

/**
 * Device kernel for clamping to the low-order viscosity
 */
kernel void entropy_visc_clamp_to_low_order_kernel(device float *reg_coeff[[ buffer(0) ]],
                                                   device const float *h[[ buffer(1) ]],
                                                   device const float *max_wave_speed[[ buffer(2) ]],
                                                   constant float &c_avisc_low[[ buffer(3) ]],
                                                   constant int &n[[ buffer(4) ]],
                                                   uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float low_order_visc = c_avisc_low * h[idx] * max_wave_speed[idx];
  reg_coeff[idx] = fmin(reg_coeff[idx], low_order_visc);
}

/**
 * Device kernel for dividing by the multiplicity (smoothing)
 */
kernel void entropy_visc_smooth_divide_kernel(device float *reg_coeff[[ buffer(0) ]],
                                              device const float *temp_field[[ buffer(1) ]],
                                              device const float *mult_field[[ buffer(2) ]],
                                              constant int &n[[ buffer(3) ]],
                                              uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  reg_coeff[idx] = temp_field[idx] / mult_field[idx];
}
