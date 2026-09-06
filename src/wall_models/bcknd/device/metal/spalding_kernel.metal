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
 * Metal compute kernel for Spalding's wall model.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Newton solver for the algebraic equation defined by Spalding's law.
 */
static float spalding_solve(const float u, const float y, const float guess,
                            const float nu, const float kappa,
                            const float B) {
  float utau = guess;
  float yp, up, f, df, old, error;
  const int maxiter = 100;

  for (int k = 0; k < maxiter; ++k) {
    up = u / utau;
    yp = y * utau / nu;
    old = utau;

    /* Evaluate function and its derivative */
    f = (up + exp(-kappa * B) *
         (exp(kappa * up) - 1.0f - kappa * up -
          0.5f * (kappa * up) * (kappa * up) -
          (1.0f / 6.0f) * (kappa * up) * (kappa * up) * (kappa * up)) -
         yp);

    df = (-y / nu - u / (utau * utau) -
          kappa * up / utau * exp(-kappa * B) *
          (exp(kappa * up) - 1.0f - kappa * up -
           0.5f * (kappa * up) * (kappa * up)));

    /* Update solution */
    utau -= f / df;

    error = fabs((old - utau) / old);

    if (error < 1e-3f) {
      break;
    }
  }

  return utau;
}

kernel void spalding_compute_kernel(device const float *u_d [[ buffer(0) ]],
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
                                    constant int &tstep [[ buffer(15) ]],
                                    uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n_nodes) return;

  const int i = (int) idx;

  float ui = u_d[i];
  float vi = v_d[i];
  float wi = w_d[i];
  const float rho = rho_w_d[i];

  /* Load normal vectors and the sampling distance once */
  const float nx = n_x_d[i];
  const float ny = n_y_d[i];
  const float nz = n_z_d[i];
  const float h = h_d[i];

  /* Project on tangential direction */
  const float normu = ui * nx + vi * ny + wi * nz;

  ui -= normu * nx;
  vi -= normu * ny;
  wi -= normu * nz;

  const float magu = sqrt(ui * ui + vi * vi + wi * wi);

  /* Get initial guess for the Newton solver */
  float guess;
  if (tstep == 1) {
    guess = sqrt(magu * nu_d[i] / h);
  } else {
    guess = tau_x_d[i] * tau_x_d[i] +
      tau_y_d[i] * tau_y_d[i] +
      tau_z_d[i] * tau_z_d[i];
    guess = sqrt(sqrt(guess));
  }

  /* Solve for utau using Newton's method */
  const float utau = spalding_solve(magu, h, guess, nu_d[i], kappa, B);

  /* Distribute according to the velocity vector */
  tau_x_d[i] = -rho * utau * utau * ui / magu;
  tau_y_d[i] = -rho * utau * utau * vi / magu;
  tau_z_d[i] = -rho * utau * utau * wi / magu;
}
