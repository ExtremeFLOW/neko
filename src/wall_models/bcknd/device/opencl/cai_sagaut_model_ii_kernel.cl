#ifndef __WALL_MODELS_CAI_SAGAUT_MODEL_II_KERNEL_CL__
#define __WALL_MODELS_CAI_SAGAUT_MODEL_II_KERNEL_CL__
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

/**
 * Device kernel for the OpenCL Cai & Sagaut Model-II wall model.
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
__kernel void cai_sagaut_model_ii_compute_kernel(
    __global const real * __restrict__ u_d,
    __global const real * __restrict__ v_d,
    __global const real * __restrict__ w_d,
    __global const real * __restrict__ n_x_d,
    __global const real * __restrict__ n_y_d,
    __global const real * __restrict__ n_z_d,
    __global const real * __restrict__ nu_d,
    __global const real * __restrict__ rho_w_d,
    __global const real * __restrict__ h_d,
    __global real * __restrict__ tau_x_d,
    __global real * __restrict__ tau_y_d,
    __global real * __restrict__ tau_z_d,
    const int n_nodes,
    const real kappa,
    const real B,
    const real p,
    const real s) {
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  const real one = (real) 1.0;
  const real half = (real) 0.5;
  const real eps = (sizeof(real) == sizeof(float)) ?
    (real) FLT_EPSILON : (real) DBL_EPSILON;
  const real e_const = exp(kappa * B);

  for (int i = idx; i < n_nodes; i += str) {
    real ui = u_d[i];
    real vi = v_d[i];
    real wi = w_d[i];
    const real rho = rho_w_d[i];
    const real nx = n_x_d[i];
    const real ny = n_y_d[i];
    const real nz = n_z_d[i];

    const real normu = ui * nx + vi * ny + wi * nz;
    ui -= normu * nx;
    vi -= normu * ny;
    wi -= normu * nz;

    const real magu = sqrt(ui * ui + vi * vi + wi * wi);

    if (magu < eps) {
      tau_x_d[i] = (real) 0.0;
      tau_y_d[i] = (real) 0.0;
      tau_z_d[i] = (real) 0.0;
      continue;
    }

    const real rey = magu * h_d[i] / nu_d[i];
    const real blend = exp(-pow(rey / s, p));
    const real warg = kappa * e_const * rey;
    const real a = one / (one + half * log(one + warg));
    real wlam = log(one + a * warg);
    wlam = wlam / (one + wlam) * (one + log(warg / wlam));

    const real up = blend * sqrt(rey) + (one - blend) * wlam / kappa;
    const real utau = magu / (up + eps);

    tau_x_d[i] = -rho * utau * utau * ui / (magu + eps);
    tau_y_d[i] = -rho * utau * utau * vi / (magu + eps);
    tau_z_d[i] = -rho * utau * utau * wi / (magu + eps);
  }
}

#endif
