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

#include <device/device_config.h>
#include <device/cuda/check.h>

#include "prs_res_kernel.h"
#include "vel_res_update_kernel.h"

extern "C" {

  void pnpn_prs_res_part1_cuda(void *ta1, void *ta2, void *ta3, 
			       void *wa1, void *wa2, void *wa3,
			       void *f_u, void *f_v, void *f_w,
			       void *B, void *h1, real *mu,
			       real *rho, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    prs_res_part1_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) ta1, (real *) ta2,
                                      (real *) ta3, (real *) wa1,
                                      (real *) wa2, (real *) wa3,
                                      (real *) f_u, (real *) f_v,
                                      (real *) f_w, (real *) B,
                                      (real *) h1, *mu, *rho, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  void pnpn_prs_res_part2_cuda(void *p_res, void *wa1, void *wa2,
			       void *wa3, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    prs_res_part2_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) p_res, (real *) wa1,
                                      (real *) wa2, (real *) wa3, *n);
    CUDA_CHECK(cudaGetLastError());
    
  }

  void pnpn_prs_res_part3_cuda(void *p_res, void *ta1, void *ta2,
			       void *ta3, real *dtbd, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    prs_res_part3_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) p_res, (real *) ta1,
                                      (real *) ta2, (real *) ta3,
                                      *dtbd, *n);
    CUDA_CHECK(cudaGetLastError());
    
  }
  
  void pnpn_vel_res_update_cuda(void *u_res, void *v_res, void *w_res,
				void *ta1, void *ta2, void *ta3,
				void *f_u, void *f_v, void *f_w, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    vel_res_update_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) u_res, (real *) v_res,
                                      (real *) w_res, (real *) ta1,
                                      (real *) ta2, (real *) ta3,
                                      (real *) f_u, (real *) f_v,
                                      (real *) f_w, *n);
    CUDA_CHECK(cudaGetLastError());
  }
    
}

