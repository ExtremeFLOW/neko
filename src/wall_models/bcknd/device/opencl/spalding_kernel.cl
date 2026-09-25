#ifndef __WALL_MODELS_SPALDING_KERNEL_CL__
#define __WALL_MODELS_SPALDING_KERNEL_CL__
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
 * Newton solver for the algebraic equation defined by Spalding's law.
 */
inline real spalding_solve(const real u, const real y, const real guess,
                           const real nu, const real kappa, const real B) {
  real utau = guess;
  real yp, up, f, df, old, error;
  const int maxiter = 100;
  const real one = (real) 1.0;

  for (int k = 0; k < maxiter; ++k) {
    up = u / utau;
    yp = y * utau / nu;
    old = utau;

    /* Evaluate function and its derivative */
    f = (up + exp(-kappa * B) *
         (exp(kappa * up) - one - kappa * up -
          (real) 0.5 * (kappa * up) * (kappa * up) -
          (one / (real) 6.0) * (kappa * up) * (kappa * up) * (kappa * up)) -
         yp);

    df = (-y / nu - u / (utau * utau) -
          kappa * up / utau * exp(-kappa * B) *
          (exp(kappa * up) - one - kappa * up -
           (real) 0.5 * (kappa * up) * (kappa * up)));

    /* Update solution */
    utau -= f / df;

    error = fabs((old - utau) / old);

    if (error < (real) 1e-3) {
      break;
    }
  }

  return utau;
}

/**
 * Device kernel for Spalding's wall model.
 */
__kernel void spalding_compute_kernel(
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
    const int tstep) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n_nodes; i += str) {
    real ui = u_d[i];
    real vi = v_d[i];
    real wi = w_d[i];
    const real rho = rho_w_d[i];

    /* Load normal vectors and the sampling distance once */
    const real nx = n_x_d[i];
    const real ny = n_y_d[i];
    const real nz = n_z_d[i];
    const real h = h_d[i];

    /* Project on tangential direction */
    const real normu = ui * nx + vi * ny + wi * nz;

    ui -= normu * nx;
    vi -= normu * ny;
    wi -= normu * nz;

    const real magu = sqrt(ui * ui + vi * vi + wi * wi);

    /* Get initial guess for the Newton solver */
    real guess;
    if (tstep == 1) {
      guess = sqrt(magu * nu_d[i] / h);
    } else {
      guess = tau_x_d[i] * tau_x_d[i] +
        tau_y_d[i] * tau_y_d[i] +
        tau_z_d[i] * tau_z_d[i];
      guess = sqrt(sqrt(guess));
    }

    /* Solve for utau using Newton's method */
    const real utau = spalding_solve(magu, h, guess, nu_d[i], kappa, B);

    /* Distribute according to the velocity vector */
    tau_x_d[i] = -rho * utau * utau * ui / magu;
    tau_y_d[i] = -rho * utau * utau * vi / magu;
    tau_z_d[i] = -rho * utau * utau * wi / magu;
  }
}

#endif // __WALL_MODELS_SPALDING_KERNEL_CL__
