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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "ale_kinematics_kernel.h"


extern "C" {
  void add_kinematics_to_mesh_velocity_hip(
      void *wx, void *wy, void *wz,
      void *x_ref, void *y_ref, void *z_ref,
      void *phi, void *x, void *y, void *z,
      kinematics_params_t kin_params, 
      int n)                              
  {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((n + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(
        HIP_KERNEL_NAME(ale_add_kinematics_kernel<real>),
        nblcks, nthrds, 0, (hipStream_t)glb_cmd_queue,
        n,
        (real*)wx, (real*)wy, (real*)wz,
        (real*)x_ref, (real*)y_ref, (real*)z_ref,
        (real*)phi, (real*)x, (real*)y, (real*)z,
        kin_params);

    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for compute_cheap_dist
   */
  void compute_cheap_dist_hip(void *d_d, void *x_d, void *y_d, void *z_d, 
            int lx, int ly, int lz, int nel, 
            int local_iters, void *nchange_d) {

    const dim3 nthrds(256, 1, 1);
    const dim3 nblcks((nel + 256 - 1) / 256, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(compute_cheap_dist_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t)glb_cmd_queue, 
                       (real*)d_d, (real*)x_d, (real*)y_d, (real*)z_d, 
                       lx, ly, lz, nel, local_iters, (int*)nchange_d);
                       
    HIP_CHECK(hipGetLastError());
  }

} // extern "C"