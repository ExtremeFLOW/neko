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
 * Metal compute kernels for the Duprat wall model.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Dimensionless velocity U*(y*) of the model of Duprat et al. (2011).
 * Integrates dU* / dy* = (s_p (1 - alpha)^(3/2) y* + alpha) / (1 + nu_t / nu),
 * Eq. (5), with the eddy viscosity of Eq. (6), using a composite 12 x 5-point
 * Gauss-Legendre rule in the stretched coordinate s = ln(1 + y*).
 * @param ys The upper integration limit y*.
 * @param alpha The ratio u_tau^2 / u_tau_p^2.
 * @param sp The sign of the pressure gradient relative to the flow.
 * @param kappa The von Karman coefficient.
 * @param beta The exponent of the pressure-gradient term.
 * @param A The damping constant.
 */
static float duprat_u_star(const float ys, const float alpha, const float sp,
                           const float kappa, const float beta, const float A) {
  const float gl_x[5] = {0.0469100770306680f, 0.2307653449471585f, 0.5f,
                         0.7692346550528415f, 0.9530899229693320f};
  const float gl_w[5] = {0.1184634425280945f, 0.2393143352496832f,
                         0.2844444444444444f, 0.2393143352496832f,
                         0.1184634425280945f};
  const int n_sub = 12;
  const float one = 1.0f;
  float us = 0.0f;

  if (ys <= 0.0f)
    return us;

  const float oma = (one - alpha > 0.0f) ? one - alpha : 0.0f;
  const float apg = pow(oma, 1.5f);
  const float damp = one + A * alpha * alpha * alpha;
  const float ds = log(one + ys) / (float)n_sub;

  for (int k = 0; k < n_sub; ++k) {
    const float s0 = (float)k * ds;
    for (int j = 0; j < 5; ++j) {
      const float eta = exp(s0 + gl_x[j] * ds) - one;
      const float edamp = one - exp(-eta / damp);
      const float nut =
          kappa * eta * pow(alpha + eta * apg, beta) * edamp * edamp;
      /* dU^+ / dy^+ times the Jacobian d(eta) / ds = 1 + eta */
      us += gl_w[j] * ds * (sp * apg * eta + alpha) / (one + nut) * (one + eta);
    }
  }
  return us;
}

/**
 * Residual f(utau) = u_tau_p U*(y*) - u of the model of Duprat et al.
 */
static float duprat_residual(const float utau, const float u, const float y,
                             const float nu, const float up, const float sp,
                             const float kappa, const float beta,
                             const float A) {
  const float utp = sqrt(utau * utau + up * up);
  /* alpha without the square root, so that it cannot exceed 1 */
  const float alpha = utau * utau / (utau * utau + up * up);
  return utp * duprat_u_star(y * utp / nu, alpha, sp, kappa, beta, A) - u;
}

/**
 * Safeguarded Newton solver for the friction velocity: the root is
 * bracketed by a sign change of the residual first, and Newton steps with a
 * finite-difference derivative fall back to bisection when they leave the
 * bracket.
 */
static float duprat_solve(const float u, const float y, const float guess,
                          const float nu, const float up, const float sp,
                          const float kappa, const float beta, const float A) {
  const int maxiter = 100;
  const float eps = FLT_EPSILON;
  const float tol = 10.0f * FLT_EPSILON;
  const float tiny = sqrt(eps) * u;

  /* Bracket the root, f(lo) < 0 < f(hi) */
  float lo = (guess > tiny) ? guess : tiny;
  float hi = lo;
  float f_lo = duprat_residual(lo, u, y, nu, up, sp, kappa, beta, A);
  float f_hi = f_lo;
  for (int k = 0; k < maxiter && f_hi <= 0.0f; ++k) {
    hi = 2.0f * hi;
    f_hi = duprat_residual(hi, u, y, nu, up, sp, kappa, beta, A);
  }
  for (int k = 0; k < maxiter && f_lo >= 0.0f && lo > tiny; ++k) {
    lo = (0.5f * lo > tiny) ? 0.5f * lo : tiny;
    f_lo = duprat_residual(lo, u, y, nu, up, sp, kappa, beta, A);
  }

  /* No sign change down to utau ~ 0: the shear stress vanishes */
  if (f_lo >= 0.0f)
    return lo;

  float utau = 0.5f * (lo + hi);
  if (guess > lo && guess < hi)
    utau = guess;
  float f = duprat_residual(utau, u, y, nu, up, sp, kappa, beta, A);

  for (int k = 0; k < maxiter; ++k) {
    /* Shrink the bracket */
    if (f < 0.0f) {
      lo = utau;
    } else {
      hi = utau;
    }

    /* Newton step, or bisection if it leaves the bracket */
    const float delta = sqrt(eps) * utau;
    const float df =
        (duprat_residual(utau + delta, u, y, nu, up, sp, kappa, beta, A) - f) /
        delta;
    float step = f / df;
    if (df <= 0.0f || utau - step <= lo || utau - step >= hi) {
      step = utau - 0.5f * (lo + hi);
    }
    utau -= step;

    if (fabs(step) < tol * utau || (hi - lo) < tol * utau)
      break;
    f = duprat_residual(utau, u, y, nu, up, sp, kappa, beta, A);
  }
  return utau;
}

/**
 * Metal kernel computing the wall shear stress with the model of Duprat et al.
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
kernel void duprat_compute_kernel(
    device const float *u_d [[buffer(0)]],
    device const float *v_d [[buffer(1)]],
    device const float *w_d [[buffer(2)]],
    device const float *n_x_d [[buffer(3)]],
    device const float *n_y_d [[buffer(4)]],
    device const float *n_z_d [[buffer(5)]],
    device const float *nu_d [[buffer(6)]],
    device const float *rho_w_d [[buffer(7)]],
    device const float *h_d [[buffer(8)]],
    device const float *dpds_d [[buffer(9)]],
    device float *tau_x_d [[buffer(10)]], device float *tau_y_d [[buffer(11)]],
    device float *tau_z_d [[buffer(12)]], device float *alpha_d [[buffer(13)]],
    constant int &n_nodes [[buffer(14)]], constant float &kappa [[buffer(15)]],
    constant float &beta [[buffer(16)]], constant float &A [[buffer(17)]],
    constant int &tstep [[buffer(18)]], uint idx [[thread_position_in_grid]]) {
  if (idx >= (uint)n_nodes)
    return;

  const int i = (int)idx;
  const float eps = FLT_EPSILON;

  float ui = u_d[i];
  float vi = v_d[i];
  float wi = w_d[i];
  const float rho = rho_w_d[i];
  const float nx = n_x_d[i];
  const float ny = n_y_d[i];
  const float nz = n_z_d[i];
  const float h = h_d[i];
  const float nu = nu_d[i];
  const float dp = dpds_d[i];

  /* Project on tangential direction */
  const float normu = ui * nx + vi * ny + wi * nz;
  ui -= normu * nx;
  vi -= normu * ny;
  wi -= normu * nz;

  const float magu = sqrt(ui * ui + vi * vi + wi * wi);

  /* No tangential velocity, no shear stress */
  if (magu <= eps) {
    tau_x_d[i] = 0.0f;
    tau_y_d[i] = 0.0f;
    tau_z_d[i] = 0.0f;
    alpha_d[i] = 1.0f;
    return;
  }

  /* Pressure velocity u_p = |(nu / rho) dp/ds|^(1/3), Eq. (4) */
  const float up = pow(nu * fabs(dp) / rho, 1.0f / 3.0f);
  const float sp = copysign(1.0f, dp);

  /* Get initial guess for the Newton solver */
  float guess = tau_x_d[i] * tau_x_d[i] + tau_y_d[i] * tau_y_d[i] +
                tau_z_d[i] * tau_z_d[i];
  if (tstep == 1 || guess <= 0.0f) {
    guess = sqrt(magu * nu / h);
  } else {
    guess = sqrt(sqrt(guess) / rho);
  }

  const float utau = duprat_solve(magu, h, guess, nu, up, sp, kappa, beta, A);

  /* Distribute according to the velocity vector */
  tau_x_d[i] = -rho * utau * utau * ui / magu;
  tau_y_d[i] = -rho * utau * utau * vi / magu;
  tau_z_d[i] = -rho * utau * utau * wi / magu;

  /* alpha = u_tau^2 / u_tau_p^2, 1 without pressure gradient */
  alpha_d[i] = (up > 0.0f) ? utau * utau / (utau * utau + up * up) : 1.0f;
}

/**
 * Metal kernel updating the filtered wall-tangential pressure gradient.
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
kernel void duprat_update_dpds_kernel(device const float *dpx_d [[buffer(0)]],
                                      device const float *dpy_d [[buffer(1)]],
                                      device const float *dpz_d [[buffer(2)]],
                                      device const float *u_d [[buffer(3)]],
                                      device const float *v_d [[buffer(4)]],
                                      device const float *w_d [[buffer(5)]],
                                      device const float *n_x_d [[buffer(6)]],
                                      device const float *n_y_d [[buffer(7)]],
                                      device const float *n_z_d [[buffer(8)]],
                                      device const float *nu_d [[buffer(9)]],
                                      device const float *rho_w_d
                                      [[buffer(10)]],
                                      device const float *h_d [[buffer(11)]],
                                      constant float &eps [[buffer(12)]],
                                      device float *dpds_d [[buffer(13)]],
                                      constant int &n_nodes [[buffer(14)]],
                                      uint idx [[thread_position_in_grid]]) {
  if (idx >= (uint)n_nodes)
    return;

  const int i = (int)idx;
  const float eps_mach = FLT_EPSILON;

  /* Wall-parallel part of the sampled velocity */
  const float normu = u_d[i] * n_x_d[i] + v_d[i] * n_y_d[i] + w_d[i] * n_z_d[i];
  const float ui = u_d[i] - normu * n_x_d[i];
  const float vi = v_d[i] - normu * n_y_d[i];
  const float wi = w_d[i] - normu * n_z_d[i];
  const float magu = sqrt(ui * ui + vi * vi + wi * wi);

  /* Without a flow direction, keep the previous filtered value */
  if (magu <= eps_mach) {
    return;
  }

  /* Streamwise pressure gradient, positive when adverse */
  float dp = (dpx_d[i] * ui + dpy_d[i] * vi + dpz_d[i] * wi) / magu;

  /* Limit to the Stokes bound, u_p <= sqrt(u nu / h) */
  const float dp_max =
      rho_w_d[i] / nu_d[i] * pow(magu * nu_d[i] / h_d[i], 1.5f);
  if (fabs(dp) > dp_max)
    dp = copysign(dp_max, dp);

  /* Exponential moving average with weight eps = 1 - exp(-dt / T) */
  dpds_d[i] = (1.0f - eps) * dpds_d[i] + eps * dp;
}
