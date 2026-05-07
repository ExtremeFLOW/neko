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

#ifndef __FLUID_ENTROPY_VISCOSITY_KERNEL_H__
#define __FLUID_ENTROPY_VISCOSITY_KERNEL_H__

#include <device/device_config.h>

/**
 * Kernel for computing entropy residual from BDF time derivative
 */
template<typename T>
__global__ void entropy_visc_compute_residual_kernel(T * __restrict__ entropy_residual,
                                                     const T * __restrict__ S,
                                                     const T * __restrict__ S_lag1,
                                                     const T * __restrict__ S_lag2,
                                                     const T * __restrict__ S_lag3,
                                                     const T bdf1,
                                                     const T bdf2,
                                                     const T bdf3,
                                                     const T bdf4,
                                                     const T dt,
                                                     const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

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
template<typename T>
__global__ void entropy_visc_compute_viscosity_kernel(T * __restrict__ reg_coeff,
                                                      const T * __restrict__ entropy_residual,
                                                      const T * __restrict__ h,
                                                      const T c_avisc_entropy,
                                                      const T n_S,
                                                      const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    reg_coeff[i] = c_avisc_entropy * h[i] * h[i] * entropy_residual[i] / n_S;
  }
}

/**
 * Kernel for applying element-wise maximum
 * Each block handles one element
 */
template<typename T>
__global__ void entropy_visc_apply_element_max_kernel(T * __restrict__ reg_coeff,
                                                      const int lx3,
                                                      const int nelv) {
  __shared__ T sdata[1024];

  const int el = blockIdx.x;
  const int tid = threadIdx.x;

  if (el >= nelv) return;

  const int offset = el * lx3;

  T max_val = -1e30;
  for (int i = tid; i < lx3; i += blockDim.x) {
    max_val = max(max_val, reg_coeff[offset + i]);
  }

  sdata[tid] = max_val;
  __syncthreads();

  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata[tid] = max(sdata[tid], sdata[tid + s]);
    }
    __syncthreads();
  }

  max_val = sdata[0];
  __syncthreads();

  for (int i = tid; i < lx3; i += blockDim.x) {
    reg_coeff[offset + i] = max_val;
  }
}

/**
 * Kernel for clamping to low-order viscosity
 */
template<typename T>
__global__ void entropy_visc_clamp_to_low_order_kernel(T * __restrict__ reg_coeff,
                                                       const T * __restrict__ h,
                                                       const T * __restrict__ max_wave_speed,
                                                       const T c_avisc_low,
                                                       const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T low_order_visc = c_avisc_low * h[i] * max_wave_speed[i];
    reg_coeff[i] = min(reg_coeff[i], low_order_visc);
  }
}

/**
 * Kernel for dividing by multiplicity (smoothing)
 */
template<typename T>
__global__ void entropy_visc_smooth_divide_kernel(T * __restrict__ reg_coeff,
                                                  const T * __restrict__ temp_field,
                                                  const T * __restrict__ mult_field,
                                                  const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    reg_coeff[i] = temp_field[i] / mult_field[i];
  }
}

#endif // __FLUID_ENTROPY_VISCOSITY_KERNEL_H__

