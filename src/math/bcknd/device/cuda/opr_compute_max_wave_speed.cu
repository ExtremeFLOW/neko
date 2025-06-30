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

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include <math/bcknd/device/device_config.h>
#include <math/bcknd/device/cuda/mathops_kernel.h>

void cuda_compute_max_wave_speed(real *max_wave_speed_d, 
                                 real *u_d, real *v_d, real *w_d,
                                 real gamma, real *p_d, real *rho_d, 
                                 int n) {
  
  const dim3 nthrds(256);
  const dim3 nblcks((n + 256 - 1) / 256);
  
#ifdef HAVE_REAL_SINGLE
  compute_max_wave_speed_kernel<float>
    <<<nblcks, nthrds, 0, (cudaStream_t) glb_cmd_strm>>>(max_wave_speed_d, 
                                                          u_d, v_d, w_d, 
                                                          gamma, p_d, rho_d, n);
#else
  compute_max_wave_speed_kernel<double>
    <<<nblcks, nthrds, 0, (cudaStream_t) glb_cmd_strm>>>(max_wave_speed_d, 
                                                          u_d, v_d, w_d, 
                                                          gamma, p_d, rho_d, n);
#endif
  CUDA_CHECK(cudaGetLastError());
  
}

#ifdef __cplusplus
}
#endif 