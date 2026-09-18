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

#include <device/device_config.h>
#include <device/cuda/check.h>

#include "cai_sagaut_model_ii_kernel.h"

extern "C" {
  /**
   * Fortran wrapper for the CUDA Cai & Sagaut Model-II kernel.
   * @param u_d The sampled x-velocity field.
   * @param v_d The sampled y-velocity field.
   * @param w_d The sampled z-velocity field.
   * @param u_d Sampled x velocity.
   * @param v_d Sampled y velocity.
   * @param w_d Sampled z velocity.
   * @param n_x_d The x-component of the wall normals.
   * @param n_y_d The y-component of the wall normals.
   * @param n_z_d The z-component of the wall normals.
   * @param nu_d The sampled kinematic viscosity at wall points.
   * @param rho_w_d The sampled density at wall points.
   * @param h_d The wall-model sampling distances.
   * @param tau_x_d The x-component of the wall shear stress.
   * @param tau_y_d The y-component of the wall shear stress.
   * @param tau_z_d The z-component of the wall shear stress.
   * @param n_nodes The number of wall points.
   * @param lx The number of GLL points per direction.
   * @param kappa The von Karman coefficient.
   * @param B The log-law intercept.
   * @param p The blending exponent.
   * @param s The blending scale.
   */
  void cuda_cai_sagaut_model_ii_compute(void *u_d, void *v_d, void *w_d,
                                       void *n_x_d, void *n_y_d, void *n_z_d,
                                       void *nu_d, void *rho_w_d, void *h_d,
                                       void *tau_x_d, void *tau_y_d,
                                       void *tau_z_d, int *n_nodes,
                                       real *kappa, real *B, real *p,
                                       real *s) {
    const dim3 nthrds(256, 1, 1);
    const dim3 nblcks(((*n_nodes) + 256 - 1) / 256, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    if (*n_nodes > 0) {
      cai_sagaut_model_ii_compute<real>
         <<<nblcks, nthrds, 0, stream>>>((real *) u_d, (real *) v_d,
                         (real *) w_d,
                         (real *) n_x_d, (real *) n_y_d,
                         (real *) n_z_d, (real *) nu_d, (real *) rho_w_d,
                         (real *) h_d, (real *) tau_x_d, (real *) tau_y_d,
                         (real *) tau_z_d, *n_nodes, *kappa, *B, *p,
                         *s);
      CUDA_CHECK(cudaGetLastError());
    }
  }
}
