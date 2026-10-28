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

#ifndef __WALL_MODELS_REICHARDT_KERNEL_H__
#define __WALL_MODELS_REICHARDT_KERNEL_H__

#include <cfloat>
#include <cmath>

/**
 * Reichardt's law of the wall, u+ as a function of y+.
 * @param yp The wall-normal distance in wall units.
 * @param kappa The von Karman coefficient.
 * @param C The amplitude of the exponential correction.
 * @param B1 The damping length scale, in wall units.
 * @param B2 The decay length scale of the second exponential term.
 */
template<typename T>
__device__ T reichardt_up(const T yp, const T kappa, const T C,
                          const T B1, const T B2) {
  const T one = static_cast<T>(1.0);

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
template<typename T>
__device__ T reichardt_dup(const T yp, const T kappa, const T C,
                           const T B1, const T B2) {
  const T one = static_cast<T>(1.0);

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
template<typename T>
__device__ T reichardt_solve(const T u, const T y, const T guess,
                             const T nu, const T kappa, const T C,
                             const T B1, const T B2, const T tol) {
  T utau = guess;
  const int maxiter = 100;

  for (int k = 0; k < maxiter; ++k) {
    const T old = utau;
    const T yp = y * utau / nu;
    const T up = reichardt_up(yp, kappa, C, B1, B2);

    // Evaluate function and its derivative
    const T f = utau * up - u;
    const T df = up + yp * reichardt_dup(yp, kappa, C, B1, B2);

    // Update solution, keeping utau positive
    utau -= f / df;
    if (utau <= static_cast<T>(0.0)) {
      utau = static_cast<T>(0.5) * old;
    }

    if (fabs((old - utau) / old) < tol) {
      break;
    }
  }

  return utau;
}

/**
 * HIP kernel for Reichardt's wall model.
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
template<typename T>
__global__ void reichardt_compute(const T * __restrict__ u_d,
                                  const T * __restrict__ v_d,
                                  const T * __restrict__ w_d,
                                  const T * __restrict__ n_x_d,
                                  const T * __restrict__ n_y_d,
                                  const T * __restrict__ n_z_d,
                                  const T * __restrict__ nu_d,
                                  const T * __restrict__ rho_w_d,
                                  const T * __restrict__ h_d,
                                  T * __restrict__ tau_x_d,
                                  T * __restrict__ tau_y_d,
                                  T * __restrict__ tau_z_d,
                                  const int n_nodes,
                                  const T kappa,
                                  const T C,
                                  const T B1,
                                  const T B2,
                                  const int tstep) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const T eps = (sizeof(T) == sizeof(float)) ?
    static_cast<T>(FLT_EPSILON) :
    static_cast<T>(DBL_EPSILON);
  const T tol = (sizeof(T) == sizeof(float)) ?
    static_cast<T>(10.0 * FLT_EPSILON) :
    static_cast<T>(1e-8);

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

    // Project on tangential direction
    const T normu = ui * nx + vi * ny + wi * nz;
    ui -= normu * nx;
    vi -= normu * ny;
    wi -= normu * nz;

    const T magu = sqrt(ui * ui + vi * vi + wi * wi);

    // No tangential velocity, no shear stress
    if (magu <= eps) {
      tau_x_d[i] = static_cast<T>(0.0);
      tau_y_d[i] = static_cast<T>(0.0);
      tau_z_d[i] = static_cast<T>(0.0);
      continue;
    }

    // Get initial guess for the Newton solver
    T guess = tau_x_d[i] * tau_x_d[i] +
      tau_y_d[i] * tau_y_d[i] +
      tau_z_d[i] * tau_z_d[i];
    if (tstep == 1 || guess <= static_cast<T>(0.0)) {
      guess = sqrt(magu * nu / h);
    } else {
      guess = sqrt(sqrt(guess) / rho);
    }

    const T utau = reichardt_solve(magu, h, guess, nu, kappa, C, B1, B2,
                                   tol);

    // Distribute according to the velocity vector
    tau_x_d[i] = -rho * utau * utau * ui / magu;
    tau_y_d[i] = -rho * utau * utau * vi / magu;
    tau_z_d[i] = -rho * utau * utau * wi / magu;
  }
}

#endif // __WALL_MODELS_REICHARDT_KERNEL_H__
