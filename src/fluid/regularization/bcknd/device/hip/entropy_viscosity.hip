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

#include <device/device_config.h>
#include <device/hip/check.h>
#include "entropy_viscosity_kernel.h"

extern "C" {

void hip_entropy_visc_compute_residual(void *entropy_residual,
                                       void *S, void *S_lag1,
                                       void *S_lag2, void *S_lag3,
                                       real *bdf1, real *bdf2,
                                       real *bdf3, real *bdf4,
                                       real *dt, int *n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  hipLaunchKernelGGL(HIP_KERNEL_NAME(entropy_visc_compute_residual_kernel<real>),
                     nblcks, nthrds, 0, stream,
                     (real *) entropy_residual,
                     (real *) S, (real *) S_lag1,
                     (real *) S_lag2, (real *) S_lag3,
                     *bdf1, *bdf2, *bdf3, *bdf4, *dt, *n);
  HIP_CHECK(hipGetLastError());
}

void hip_entropy_visc_compute_viscosity(void *reg_coeff,
                                        void *entropy_residual,
                                        void *h,
                                        real *c_entropy,
                                        real *n_S,
                                        int *n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  hipLaunchKernelGGL(HIP_KERNEL_NAME(entropy_visc_compute_viscosity_kernel<real>),
                     nblcks, nthrds, 0, stream,
                     (real *) reg_coeff,
                     (real *) entropy_residual,
                     (real *) h,
                     *c_entropy, *n_S, *n);
  HIP_CHECK(hipGetLastError());
}

void hip_entropy_visc_apply_element_max(void *reg_coeff,
                                        int *lx,
                                        int *nelv) {
  const int lx3 = (*lx) * (*lx) * (*lx);
  const int nthrds_x = min(1024, lx3);
  const dim3 nthrds(nthrds_x, 1, 1);
  const dim3 nblcks(*nelv, 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  hipLaunchKernelGGL(HIP_KERNEL_NAME(entropy_visc_apply_element_max_kernel<real>),
                     nblcks, nthrds, 0, stream,
                     (real *) reg_coeff, lx3, *nelv);
  HIP_CHECK(hipGetLastError());
}

void hip_entropy_visc_clamp_to_low_order(void *reg_coeff,
                                         void *h,
                                         void *max_wave_speed,
                                         real *c_max,
                                         int *n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  hipLaunchKernelGGL(HIP_KERNEL_NAME(entropy_visc_clamp_to_low_order_kernel<real>),
                     nblcks, nthrds, 0, stream,
                     (real *) reg_coeff,
                     (real *) h,
                     (real *) max_wave_speed,
                     *c_max, *n);
  HIP_CHECK(hipGetLastError());
}

void hip_entropy_visc_smooth_divide(void *reg_coeff,
                                    void *temp_field,
                                    void *mult_field,
                                    int *n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  hipLaunchKernelGGL(HIP_KERNEL_NAME(entropy_visc_smooth_divide_kernel<real>),
                     nblcks, nthrds, 0, stream,
                     (real *) reg_coeff,
                     (real *) temp_field,
                     (real *) mult_field,
                     *n);
  HIP_CHECK(hipGetLastError());
}

}

