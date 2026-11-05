#ifndef __WALL_MODELS_REICHARDT_KERNEL_CL__
#define __WALL_MODELS_REICHARDT_KERNEL_CL__
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
 * Reichardt's law of the wall, u+ as a function of y+.
 * @param yp The wall-normal distance in wall units.
 * @param kappa The von Karman coefficient.
 * @param C The amplitude of the exponential correction.
 * @param B1 The damping length scale, in wall units.
 * @param B2 The decay length scale of the second exponential term.
 */
inline real reichardt_up(const real yp, const real kappa,
                         const real C, const real B1,
                         const real B2) {
  const real one = (real) 1.0;

  return log(one + kappa * yp) / kappa +
    C * (one - exp(-yp / B1) - yp / B1 * exp(-yp / B2));
}

/**
 * Derivative du+/dy+ of Reichardt's law of the wall.
 * @param yp The wall-normal distance in wall units.
 * @param kappa The von Karman coefficient.
 * @param C The amplitude of the exponential correction.
 * @param B1 The damping length scale, in wall units.
 * @param B2 The decay length scale of the second exponential term.
 */
inline real reichardt_dup(const real yp, const real kappa,
                          const real C, const real B1,
                          const real B2) {
  const real one = (real) 1.0;

  return one / (one + kappa * yp) +
    C / B1 * (exp(-yp / B1) - (one - yp / B2) * exp(-yp / B2));
}

/**
 * Newton solver for f(utau) = utau * u+(y utau / nu) - u = 0.
 * @param u The tangential velocity magnitude.
 * @param y The wall-normal distance.
 * @param guess Initial guess.
 * @param nu The kinematic viscosity.
 * @param kappa The von Karman coefficient.
 * @param C The amplitude of the exponential correction.
 * @param B1 The damping length scale, in wall units.
 * @param B2 The decay length scale of the second exponential term.
 * @param tol The relative convergence tolerance.
 */
inline real reichardt_solve(const real u, const real y,
                            const real guess, const real nu,
                            const real kappa, const real C,
                            const real B1, const real B2,
                            const real tol) {
  real utau = guess;
  const int maxiter = 100;

  for (int k = 0; k < maxiter; ++k) {
    const real old = utau;
    const real yp = y * utau / nu;
    const real up = reichardt_up(yp, kappa, C, B1, B2);

    /* Evaluate function and its derivative */
    const real f = utau * up - u;
    const real df = up + yp * reichardt_dup(yp, kappa, C, B1, B2);

    /* Update solution, keeping utau positive */
    utau -= f / df;
    if (utau <= (real) 0.0) {
      utau = (real) 0.5 * old;
    }

    if (fabs((old - utau) / old) < tol) {
      break;
    }
  }

  return utau;
}

/**
 * OpenCL kernel for Reichardt's wall model.
 * @param u_d The sampled x-velocity.
 * @param v_d The sampled y-velocity.
 * @param w_d The sampled z-velocity.
 * @param n_x_d The x-component of the wall normals.
 * @param n_y_d The y-component of the wall normals.
 * @param n_z_d The z-component of the wall normals.
 * @param nu_d The kinematic viscosity at wall points.
 * @param rho_w_d The density at wall points.
 * @param h_d The wall-model sampling distances.
 * @param tau_x_d The x-component of the wall shear stress.
 * @param tau_y_d The y-component of the wall shear stress.
 * @param tau_z_d The z-component of the wall shear stress.
 * @param n_nodes The number of wall points.
 * @param kappa The von Karman coefficient.
 * @param C The amplitude of the exponential correction.
 * @param B1 The damping length scale, in wall units.
 * @param B2 The decay length scale of the second exponential term.
 * @param tstep The current time-step.
 */
__kernel void reichardt_compute_kernel(
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
    const real C,
    const real B1,
    const real B2,
    const int tstep) {
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  const real eps = (sizeof(real) == sizeof(float)) ?
    (real) FLT_EPSILON : (real) DBL_EPSILON;
  const real tol = (sizeof(real) == sizeof(float)) ?
    (real) (10.0f * FLT_EPSILON) : (real) 1e-8;

  for (int i = idx; i < n_nodes; i += str) {
    real ui = u_d[i];
    real vi = v_d[i];
    real wi = w_d[i];
    const real rho = rho_w_d[i];
    const real nx = n_x_d[i];
    const real ny = n_y_d[i];
    const real nz = n_z_d[i];
    const real h = h_d[i];
    const real nu = nu_d[i];

    /* Project on tangential direction */
    const real normu = ui * nx + vi * ny + wi * nz;
    ui -= normu * nx;
    vi -= normu * ny;
    wi -= normu * nz;

    const real magu = sqrt(ui * ui + vi * vi + wi * wi);

    /* No tangential velocity, no shear stress */
    if (magu <= eps) {
      tau_x_d[i] = (real) 0.0;
      tau_y_d[i] = (real) 0.0;
      tau_z_d[i] = (real) 0.0;
      continue;
    }

    /* Get initial guess for the Newton solver */
    real guess = tau_x_d[i] * tau_x_d[i] +
      tau_y_d[i] * tau_y_d[i] +
      tau_z_d[i] * tau_z_d[i];
    if (tstep == 1 || guess <= (real) 0.0) {
      guess = sqrt(magu * nu / h);
    } else {
      guess = sqrt(sqrt(guess) / rho);
    }

    const real utau = reichardt_solve(magu, h, guess, nu, kappa, C, B1, B2,
                                   tol);

    /* Distribute according to the velocity vector */
    tau_x_d[i] = -rho * utau * utau * ui / magu;
    tau_y_d[i] = -rho * utau * utau * vi / magu;
    tau_z_d[i] = -rho * utau * utau * wi / magu;
  }
}

#endif // __WALL_MODELS_REICHARDT_KERNEL_CL__
