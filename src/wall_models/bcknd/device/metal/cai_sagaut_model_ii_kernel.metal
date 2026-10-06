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
 * Metal compute kernel for the Cai & Sagaut Model-II wall model.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

kernel void cai_sagaut_model_ii_compute_kernel(device const float *u_d [[ buffer(0) ]],
                                               device const float *v_d [[ buffer(1) ]],
                                               device const float *w_d [[ buffer(2) ]],
                                               device const float *n_x_d [[ buffer(3) ]],
                                               device const float *n_y_d [[ buffer(4) ]],
                                               device const float *n_z_d [[ buffer(5) ]],
                                               device const float *nu_d [[ buffer(6) ]],
                                               device const float *rho_w_d [[ buffer(7) ]],
                                               device const float *h_d [[ buffer(8) ]],
                                               device float *tau_x_d [[ buffer(9) ]],
                                               device float *tau_y_d [[ buffer(10) ]],
                                               device float *tau_z_d [[ buffer(11) ]],
                                               constant int &n_nodes [[ buffer(12) ]],
                                               constant float &kappa [[ buffer(13) ]],
                                               constant float &B [[ buffer(14) ]],
                                               constant float &p [[ buffer(15) ]],
                                               constant float &s [[ buffer(16) ]],
                                               uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n_nodes) return;

  const int i = (int) idx;

  const float eps = FLT_EPSILON;
  const float e_const = exp(kappa * B);

  float ui = u_d[i];
  float vi = v_d[i];
  float wi = w_d[i];
  const float rho = rho_w_d[i];
  const float nx = n_x_d[i];
  const float ny = n_y_d[i];
  const float nz = n_z_d[i];

  const float normu = ui * nx + vi * ny + wi * nz;
  ui -= normu * nx;
  vi -= normu * ny;
  wi -= normu * nz;

  const float magu = sqrt(ui * ui + vi * vi + wi * wi);

  if (magu < eps) {
    tau_x_d[i] = 0.0f;
    tau_y_d[i] = 0.0f;
    tau_z_d[i] = 0.0f;
    return;
  }

  const float rey = magu * h_d[i] / nu_d[i];
  const float blend = exp(-pow(rey / s, p));
  const float warg = kappa * e_const * rey;
  const float a = 1.0f / (1.0f + 0.5f * log(1.0f + warg));
  float wlam = log(1.0f + a * warg);
  wlam = wlam / (1.0f + wlam) * (1.0f + log(warg / wlam));

  const float up = blend * sqrt(rey) + (1.0f - blend) * wlam / kappa;
  const float utau = magu / (up + eps);

  tau_x_d[i] = -rho * utau * utau * ui / (magu + eps);
  tau_y_d[i] = -rho * utau * utau * vi / (magu + eps);
  tau_z_d[i] = -rho * utau * utau * wi / (magu + eps);
}
