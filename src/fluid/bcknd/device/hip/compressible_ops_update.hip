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

#include <stdio.h>
#include <stdlib.h>
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include <fluid/bcknd/device/hip/compressible_ops_kernel.h>

extern "C" {

void hip_update_uvw(void *u, void *v, void *w,
                     void *m_x, void *m_y, void *m_z,
                     void *rho, int *n) {
  
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;
  
  hipLaunchKernelGGL(HIP_KERNEL_NAME(update_uvw_kernel<real>),
                     nblcks, nthrds, 0, stream,
                     (real *) u, (real *) v, (real *) w,
                     (real *) m_x, (real *) m_y, (real *) m_z,
                     (real *) rho, *n);

  HIP_CHECK(hipGetLastError());
  
}

void hip_update_mxyz_p_ruvw(void *m_x, void *m_y, void *m_z,
                             void *p, void *ruvw, 
                             void *u, void *v, void *w, void *E,                     
                             void *rho, real *gamma, int *n) {
  
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;
  
  hipLaunchKernelGGL(HIP_KERNEL_NAME(update_mxyz_p_ruvw_kernel<real>),
                     nblcks, nthrds, 0, stream,
                     (real *) m_x, (real *) m_y, (real *) m_z,
                     (real *) p, (real *) ruvw,
                     (real *) u, (real *) v, (real *) w,
                     (real *) E, (real *) rho, *gamma, *n);

  HIP_CHECK(hipGetLastError());
  
}

void hip_update_e(void *E, void *p, void *ruvw,
                   real *gamma, int *n) {
  
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;
  
  hipLaunchKernelGGL(HIP_KERNEL_NAME(update_e_kernel<real>),
                    nblcks, nthrds, 0, stream,
                    (real *) E, (real *) p, (real *) ruvw,
                    *gamma, *n);

  HIP_CHECK(hipGetLastError());
  
}

} 
