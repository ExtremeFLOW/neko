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

/**
 * Kernel for computing entropy residual from BDF time derivative
 */
__kernel void entropy_visc_compute_residual_kernel(__global real * __restrict__ entropy_residual,
                                                   __global const real * __restrict__ S,
                                                   __global const real * __restrict__ S_lag1,
                                                   __global const real * __restrict__ S_lag2,
                                                   __global const real * __restrict__ S_lag3,
                                                   const real bdf1,
                                                   const real bdf2,
                                                   const real bdf3,
                                                   const real bdf4,
                                                   const real dt,
                                                   const int n) {
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    entropy_residual[i] = (bdf1 * S[i]
                           - bdf2 * S_lag1[i]
                           - bdf3 * S_lag2[i]
                           - bdf4 * S_lag3[i]) / dt;
  }
}

/**
 * Kernel for computing viscosity from entropy residual
 */
__kernel void entropy_visc_compute_viscosity_kernel(__global real * __restrict__ reg_coeff,
                                                    __global const real * __restrict__ entropy_residual,
                                                    __global const real * __restrict__ h,
                                                    const real c_entropy,
                                                    const real n_S,
                                                    const int n) {
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    reg_coeff[i] = c_entropy * h[i] * h[i] * entropy_residual[i] / n_S;
  }
}

/**
 * Kernel for applying element-wise maximum
 * Each work-group handles one element
 */
__kernel void entropy_visc_apply_element_max_kernel(__global real * __restrict__ reg_coeff,
                                                    const int lx3,
                                                    const int nelv) {
  __local real sdata[1024];

  const int el = get_group_id(0);
  const int tid = get_local_id(0);
  const int local_size = get_local_size(0);

  if (el >= nelv) return;

  const int offset = el * lx3;

  real max_val = -1e30;
  for (int i = tid; i < lx3; i += local_size) {
    max_val = fmax(max_val, reg_coeff[offset + i]);
  }

  sdata[tid] = max_val;
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int s = local_size / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  max_val = sdata[0];
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int i = tid; i < lx3; i += local_size) {
    reg_coeff[offset + i] = max_val;
  }
}

/**
 * Kernel for clamping to low-order viscosity
 */
__kernel void entropy_visc_clamp_to_low_order_kernel(__global real * __restrict__ reg_coeff,
                                                     __global const real * __restrict__ h,
                                                     __global const real * __restrict__ max_wave_speed,
                                                     const real c_max,
                                                     const int n) {
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    real low_order_visc = c_max * h[i] * max_wave_speed[i];
    reg_coeff[i] = fmin(reg_coeff[i], low_order_visc);
  }
}

/**
 * Kernel for setting low-order viscosity directly (no clamping)
 */
__kernel void entropy_visc_set_low_order_kernel(__global real * __restrict__ reg_coeff,
                                                __global const real * __restrict__ h,
                                                __global const real * __restrict__ max_wave_speed,
                                                const real c_max,
                                                const int n) {
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    reg_coeff[i] = c_max * h[i] * max_wave_speed[i];
  }
}

/**
 * Kernel for dividing by multiplicity (smoothing)
 */
__kernel void entropy_visc_smooth_divide_kernel(__global real * __restrict__ reg_coeff,
                                                __global const real * __restrict__ temp_field,
                                                __global const real * __restrict__ mult_field,
                                                const int n) {
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    reg_coeff[i] = temp_field[i] / mult_field[i];
  }
}

