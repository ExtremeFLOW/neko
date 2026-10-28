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

#ifndef __WALL_MODELS_DUPRAT_KERNEL_H__
#define __WALL_MODELS_DUPRAT_KERNEL_H__

#include <cfloat>
#include <cmath>

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
template <typename T>
__device__ T duprat_u_star(const T ys, const T alpha, const T sp, const T kappa,
                           const T beta, const T A) {
  const T gl_x[5] = {static_cast<T>(0.0469100770306680),
                     static_cast<T>(0.2307653449471585), static_cast<T>(0.5),
                     static_cast<T>(0.7692346550528415),
                     static_cast<T>(0.9530899229693320)};
  const T gl_w[5] = {
      static_cast<T>(0.1184634425280945), static_cast<T>(0.2393143352496832),
      static_cast<T>(0.2844444444444444), static_cast<T>(0.2393143352496832),
      static_cast<T>(0.1184634425280945)};
  const int n_sub = 12;
  const T one = static_cast<T>(1.0);
  T us = static_cast<T>(0.0);

  if (ys <= static_cast<T>(0.0))
    return us;

  const T oma =
      (one - alpha > static_cast<T>(0.0)) ? one - alpha : static_cast<T>(0.0);
  const T apg = pow(oma, static_cast<T>(1.5));
  const T damp = one + A * alpha * alpha * alpha;
  const T ds = log(one + ys) / (T)n_sub;

  for (int k = 0; k < n_sub; ++k) {
    const T s0 = (T)k * ds;
    for (int j = 0; j < 5; ++j) {
      const T eta = exp(s0 + gl_x[j] * ds) - one;
      const T edamp = one - exp(-eta / damp);
      const T nut = kappa * eta * pow(alpha + eta * apg, beta) * edamp * edamp;
      /* dU^+ / dy^+ times the Jacobian d(eta) / ds = 1 + eta */
      us += gl_w[j] * ds * (sp * apg * eta + alpha) / (one + nut) * (one + eta);
    }
  }
  return us;
}

/**
 * Residual f(utau) = u_tau_p U*(y*) - u of the model of Duprat et al.
 */
template <typename T>
__device__ T duprat_residual(const T utau, const T u, const T y, const T nu,
                             const T up, const T sp, const T kappa,
                             const T beta, const T A) {
  const T utp = sqrt(utau * utau + up * up);
  /* alpha without the square root, so that it cannot exceed 1 */
  const T alpha = utau * utau / (utau * utau + up * up);
  return utp * duprat_u_star(y * utp / nu, alpha, sp, kappa, beta, A) - u;
}

/**
 * Safeguarded Newton solver for the friction velocity: the root is
 * bracketed by a sign change of the residual first, and Newton steps with a
 * finite-difference derivative fall back to bisection when they leave the
 * bracket.
 */
template <typename T>
__device__ T duprat_solve(const T u, const T y, const T guess, const T nu,
                          const T up, const T sp, const T kappa, const T beta,
                          const T A) {
  const int maxiter = 100;
  const T eps = (sizeof(T) == sizeof(float)) ? static_cast<T>(FLT_EPSILON)
                                             : static_cast<T>(DBL_EPSILON);
  const T tol = (sizeof(T) == sizeof(float))
                    ? static_cast<T>(10.0 * FLT_EPSILON)
                    : static_cast<T>(1e-8);
  const T tiny = sqrt(eps) * u;

  /* Bracket the root, f(lo) < 0 < f(hi) */
  T lo = (guess > tiny) ? guess : tiny;
  T hi = lo;
  T f_lo = duprat_residual(lo, u, y, nu, up, sp, kappa, beta, A);
  T f_hi = f_lo;
  for (int k = 0; k < maxiter && f_hi <= static_cast<T>(0.0); ++k) {
    hi = static_cast<T>(2.0) * hi;
    f_hi = duprat_residual(hi, u, y, nu, up, sp, kappa, beta, A);
  }
  for (int k = 0; k < maxiter && f_lo >= static_cast<T>(0.0) && lo > tiny;
       ++k) {
    lo = (static_cast<T>(0.5) * lo > tiny) ? static_cast<T>(0.5) * lo : tiny;
    f_lo = duprat_residual(lo, u, y, nu, up, sp, kappa, beta, A);
  }

  /* No sign change down to utau ~ 0: the shear stress vanishes */
  if (f_lo >= static_cast<T>(0.0))
    return lo;

  T utau = static_cast<T>(0.5) * (lo + hi);
  if (guess > lo && guess < hi)
    utau = guess;
  T f = duprat_residual(utau, u, y, nu, up, sp, kappa, beta, A);

  for (int k = 0; k < maxiter; ++k) {
    /* Shrink the bracket */
    if (f < static_cast<T>(0.0)) {
      lo = utau;
    } else {
      hi = utau;
    }

    /* Newton step, or bisection if it leaves the bracket */
    const T delta = sqrt(eps) * utau;
    const T df =
        (duprat_residual(utau + delta, u, y, nu, up, sp, kappa, beta, A) - f) /
        delta;
    T step = f / df;
    if (df <= static_cast<T>(0.0) || utau - step <= lo || utau - step >= hi) {
      step = utau - static_cast<T>(0.5) * (lo + hi);
    }
    utau -= step;

    if (fabs(step) < tol * utau || (hi - lo) < tol * utau)
      break;
    f = duprat_residual(utau, u, y, nu, up, sp, kappa, beta, A);
  }
  return utau;
}

/**
 * CUDA kernel computing the wall shear stress with the model of Duprat et al.
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
template <typename T>
__global__ void
duprat_compute(const T *__restrict__ u_d, const T *__restrict__ v_d,
               const T *__restrict__ w_d, const T *__restrict__ n_x_d,
               const T *__restrict__ n_y_d, const T *__restrict__ n_z_d,
               const T *__restrict__ nu_d, const T *__restrict__ rho_w_d,
               const T *__restrict__ h_d, const T *__restrict__ dpds_d,
               T *__restrict__ tau_x_d, T *__restrict__ tau_y_d,
               T *__restrict__ tau_z_d, T *__restrict__ alpha_d,
               const int n_nodes, const T kappa, const T beta, const T A,
               const int tstep) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const T eps = (sizeof(T) == sizeof(float)) ? static_cast<T>(FLT_EPSILON)
                                             : static_cast<T>(DBL_EPSILON);

  for (int i = idx; i < n_nodes; i += str) {
    T ui = u_d[i];
    T vi = v_d[i];
    T wi = w_d[i];
    const T rho = rho_w_d[i];
    const T nx = n_x_d[i];
    const T ny = n_y_d[i];
    const T nz = n_z_d[i];
    const T h = h_d[i];
    const T nu = nu_d[i];
    const T dp = dpds_d[i];

    /* Project on tangential direction */
    const T normu = ui * nx + vi * ny + wi * nz;
    ui -= normu * nx;
    vi -= normu * ny;
    wi -= normu * nz;

    const T magu = sqrt(ui * ui + vi * vi + wi * wi);

    /* No tangential velocity, no shear stress */
    if (magu <= eps) {
      tau_x_d[i] = static_cast<T>(0.0);
      tau_y_d[i] = static_cast<T>(0.0);
      tau_z_d[i] = static_cast<T>(0.0);
      alpha_d[i] = static_cast<T>(1.0);
      continue;
    }

    /* Pressure velocity u_p = |(nu / rho) dp/ds|^(1/3), Eq. (4) */
    const T up =
        pow(nu * fabs(dp) / rho, static_cast<T>(1.0) / static_cast<T>(3.0));
    const T sp = copysign(static_cast<T>(1.0), dp);

    /* Get initial guess for the Newton solver */
    T guess = tau_x_d[i] * tau_x_d[i] + tau_y_d[i] * tau_y_d[i] +
              tau_z_d[i] * tau_z_d[i];
    if (tstep == 1 || guess <= static_cast<T>(0.0)) {
      guess = sqrt(magu * nu / h);
    } else {
      guess = sqrt(sqrt(guess) / rho);
    }

    const T utau = duprat_solve(magu, h, guess, nu, up, sp, kappa, beta, A);

    /* Distribute according to the velocity vector */
    tau_x_d[i] = -rho * utau * utau * ui / magu;
    tau_y_d[i] = -rho * utau * utau * vi / magu;
    tau_z_d[i] = -rho * utau * utau * wi / magu;

    /* alpha = u_tau^2 / u_tau_p^2, 1 without pressure gradient */
    alpha_d[i] = (up > static_cast<T>(0.0))
                     ? utau * utau / (utau * utau + up * up)
                     : static_cast<T>(1.0);
  }
}

/**
 * CUDA kernel updating the filtered wall-tangential pressure gradient:
 * projection on the flow direction, Stokes limit and exponential moving
 * average.
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
template <typename T>
__global__ void
duprat_update_dpds(const T *__restrict__ dpx_d, const T *__restrict__ dpy_d,
                   const T *__restrict__ dpz_d, const T *__restrict__ u_d,
                   const T *__restrict__ v_d, const T *__restrict__ w_d,
                   const T *__restrict__ n_x_d, const T *__restrict__ n_y_d,
                   const T *__restrict__ n_z_d, const T *__restrict__ nu_d,
                   const T *__restrict__ rho_w_d, const T *__restrict__ h_d,
                   const T eps, T *__restrict__ dpds_d, const int n_nodes) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const T eps_mach = (sizeof(T) == sizeof(float)) ? static_cast<T>(FLT_EPSILON)
                                                  : static_cast<T>(DBL_EPSILON);

  for (int i = idx; i < n_nodes; i += str) {
    /* Wall-parallel part of the sampled velocity */
    const T normu = u_d[i] * n_x_d[i] + v_d[i] * n_y_d[i] + w_d[i] * n_z_d[i];
    const T ui = u_d[i] - normu * n_x_d[i];
    const T vi = v_d[i] - normu * n_y_d[i];
    const T wi = w_d[i] - normu * n_z_d[i];
    const T magu = sqrt(ui * ui + vi * vi + wi * wi);

    /* Without a flow direction, keep the previous filtered value */
    if (magu <= eps_mach) {
      continue;
    }

    /* Streamwise pressure gradient, positive when adverse */
    T dp = (dpx_d[i] * ui + dpy_d[i] * vi + dpz_d[i] * wi) / magu;

    /* Limit to the Stokes bound, u_p <= sqrt(u nu / h) */
    const T dp_max = rho_w_d[i] / nu_d[i] *
                     pow(magu * nu_d[i] / h_d[i], static_cast<T>(1.5));
    if (fabs(dp) > dp_max)
      dp = copysign(dp_max, dp);

    /* Exponential moving average with weight eps = 1 - exp(-dt / T) */
    dpds_d[i] = (static_cast<T>(1.0) - eps) * dpds_d[i] + eps * dp;
  }
}

#endif // __WALL_MODELS_DUPRAT_KERNEL_H__
