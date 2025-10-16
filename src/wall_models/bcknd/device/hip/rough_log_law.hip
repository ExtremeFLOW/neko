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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "rough_log_law_kernel.h"

extern "C" {
  void hip_rough_log_law_compute(void *u_d, void *v_d, void *w_d,
          void *ind_r_d, void *ind_s_d, void *ind_t_d, void *ind_e_d,
          void *n_x_d, void *n_y_d, void *n_z_d, void *h_d,
          void *tau_x_d, void *tau_y_d, void *tau_z_d,
          int *n_nodes, int *lx, real *kappa, real *B, real *z0, int *tstep) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n_nodes) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    rough_log_law_compute<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) u_d, (real *) v_d, (real *) w_d,
                                    (int *) ind_r_d, (int *) ind_s_d, 
                                    (int *) ind_t_d, (int *) ind_e_d,
                                    (real *) n_x_d, (real *) n_y_d, 
                                    (real *) n_z_d, (real *) h_d,
                                    (real *) tau_x_d, (real *) tau_y_d, 
                                    (real *) tau_z_d,
                                    *n_nodes, *lx, *kappa, *B, *z0);
    HIP_CHECK(hipGetLastError());
  }
}