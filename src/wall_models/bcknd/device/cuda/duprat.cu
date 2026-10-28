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

#include "duprat_kernel.h"

extern "C" {
/**
 * Fortran wrapper for the CUDA kernel of the Duprat wall model.
 * @param u_d The x component of the sampled velocity.
 * @param v_d The y component of the sampled velocity.
 * @param w_d The z component of the sampled velocity.
 * @param n_x_d The x component of the wall normals.
 * @param n_y_d The y component of the wall normals.
 * @param n_z_d The z component of the wall normals.
 * @param nu_d The kinematic viscosity at wall points.
 * @param rho_w_d The density at wall points.
 * @param h_d The wall-model sampling distances.
 * @param dpds_d The wall-tangential pressure gradient (positive if adverse).
 * @param tau_x_d The x component of the wall shear stress.
 * @param tau_y_d The y component of the wall shear stress.
 * @param tau_z_d The z component of the wall shear stress.
 * @param alpha_d The ratio u_tau^2 / u_tau_p^2, for diagnostics.
 * @param n_nodes The number of wall points.
 * @param kappa The von Karman coefficient.
 * @param beta The exponent of the pressure-gradient term.
 * @param A The damping constant.
 * @param tstep The current time-step.
 */
void cuda_duprat_compute(void *u_d, void *v_d, void *w_d, void *n_x_d,
                         void *n_y_d, void *n_z_d, void *nu_d, void *rho_w_d,
                         void *h_d, void *dpds_d, void *tau_x_d, void *tau_y_d,
                         void *tau_z_d, void *alpha_d, int *n_nodes,
                         real *kappa, real *beta, real *A, int *tstep) {
  const dim3 nthrds(256, 1, 1);
  const dim3 nblcks(((*n_nodes) + 256 - 1) / 256, 1, 1);
  const cudaStream_t stream = (cudaStream_t)glb_cmd_queue;

  if (*n_nodes > 0) {
    duprat_compute<real><<<nblcks, nthrds, 0, stream>>>(
        (real *)u_d, (real *)v_d, (real *)w_d, (real *)n_x_d, (real *)n_y_d,
        (real *)n_z_d, (real *)nu_d, (real *)rho_w_d, (real *)h_d,
        (real *)dpds_d, (real *)tau_x_d, (real *)tau_y_d, (real *)tau_z_d,
        (real *)alpha_d, *n_nodes, *kappa, *beta, *A, *tstep);
    CUDA_CHECK(cudaGetLastError());
  }
}

/**
 * Fortran wrapper for the CUDA pressure-gradient update of the Duprat
 * wall model.
 * @param dpx_d The x component of the sampled pressure gradient.
 * @param dpy_d The y component of the sampled pressure gradient.
 * @param dpz_d The z component of the sampled pressure gradient.
 * @param u_d The x component of the sampled velocity.
 * @param v_d The y component of the sampled velocity.
 * @param w_d The z component of the sampled velocity.
 * @param n_x_d The x component of the wall normals.
 * @param n_y_d The y component of the wall normals.
 * @param n_z_d The z component of the wall normals.
 * @param nu_d The kinematic viscosity at wall points.
 * @param rho_w_d The density at wall points.
 * @param h_d The wall-model sampling distances.
 * @param eps The filter weight of the current gradient.
 * @param dpds_d The filtered pressure gradient, updated in place.
 * @param n_nodes The number of wall points.
 */
void cuda_duprat_update_dpds(void *dpx_d, void *dpy_d, void *dpz_d, void *u_d,
                             void *v_d, void *w_d, void *n_x_d, void *n_y_d,
                             void *n_z_d, void *nu_d, void *rho_w_d, void *h_d,
                             real *eps, void *dpds_d, int *n_nodes) {
  const dim3 nthrds(256, 1, 1);
  const dim3 nblcks(((*n_nodes) + 256 - 1) / 256, 1, 1);
  const cudaStream_t stream = (cudaStream_t)glb_cmd_queue;

  if (*n_nodes > 0) {
    duprat_update_dpds<real><<<nblcks, nthrds, 0, stream>>>(
        (real *)dpx_d, (real *)dpy_d, (real *)dpz_d, (real *)u_d, (real *)v_d,
        (real *)w_d, (real *)n_x_d, (real *)n_y_d, (real *)n_z_d, (real *)nu_d,
        (real *)rho_w_d, (real *)h_d, *eps, (real *)dpds_d, *n_nodes);
    CUDA_CHECK(cudaGetLastError());
  }
}
}
