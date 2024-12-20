/*
 Copyright (c) 2022, The Neko Authors
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

#include "euler_rhs_kernel.h"

extern "C" {

  void euler_res_part_visc_hip(void *rhs_rho, void *Binv,
                           void *lap_rho, real *c_avisc, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME( euler_res_part_visc_kernel<real> ),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) rhs_rho, (real *) Binv, 
                       (real *) lap_rho, *c_avisc, *n);
    HIP_CHECK(hipGetLastError());
  }

  void euler_res_part_mx_flux_hip(void *f_x, void *f_y, void *f_z,
                           void *m_x, void *m_y, void *m_z,
                           void *rho_field, void *p, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME( euler_res_part_mx_flux_kernel<real> ),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) f_x, (real *) f_y, (real *) f_z,
                       (real *) m_x, (real *) m_y, (real *) m_z,
                       (real *) rho_field, (real *) p, *n);
    HIP_CHECK(hipGetLastError());
  }

  void euler_res_part_my_flux_hip(void *f_x, void *f_y, void *f_z,
                           void *m_x, void *m_y, void *m_z,
                           void *rho_field, void *p, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME( euler_res_part_my_flux_kernel<real> ),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) f_x, (real *) f_y, (real *) f_z,
                       (real *) m_x, (real *) m_y, (real *) m_z,
                       (real *) rho_field, (real *) p, *n);
    HIP_CHECK(hipGetLastError());
  }

  void euler_res_part_mz_flux_hip(void *f_x, void *f_y, void *f_z,
                           void *m_x, void *m_y, void *m_z,
                           void *rho_field, void *p, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME( euler_res_part_mz_flux_kernel<real> ),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) f_x, (real *) f_y, (real *) f_z,
                       (real *) m_x, (real *) m_y, (real *) m_z,
                       (real *) rho_field, (real *) p, *n);
    HIP_CHECK(hipGetLastError());
  }

  void euler_res_part_E_flux_hip(void *f_x, void *f_y, void *f_z,
                           void *m_x, void *m_y, void *m_z,
                           void *rho_field, void *p, void *E,
                           int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME( euler_res_part_E_flux_kernel<real> ),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) f_x, (real *) f_y, (real *) f_z,
                       (real *) m_x, (real *) m_y, (real *) m_z,
                       (real *) rho_field, (real *) p, (real *) E, *n);
    HIP_CHECK(hipGetLastError());
  }
}

