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

void hip_compute_max_wave_speed(void *max_wave_speed_d, 
                                void *u_d, void *v_d, void *w_d,
                                real *gamma, void *p_d, void *rho_d, 
                                int *n) {
  
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  
  hipLaunchKernelGGL(compute_max_wave_speed_kernel<real>, nblcks, nthrds, 0,
                     (hipStream_t) glb_cmd_queue, (real *) max_wave_speed_d, 
                     (real *) u_d, (real *) v_d, (real *) w_d, *gamma, (real *) p_d, (real *) rho_d, *n);
  HIP_CHECK(hipGetLastError());
  
}

} 