#ifndef __WALL_MODELS_RICHARDSON_KERNEL_CL__
#define __WALL_MODELS_RICHARDSON_KERNEL_CL__
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

/*
 * Similarity laws and corrections for the STABLE regime:
 * Based on Mauritsen et al. 2007
 */

inline real f_tau_stable(const real Ri_b) {
  return (real) 0.17 * ((real) 0.25 +
                        (real) 0.75 / ((real) 1.0 + (real) 4.0 * Ri_b));
}

inline real f_theta_stable(const real Ri_b) {
  return (real) -0.145 / ((real) 1.0 + (real) 4.0 * Ri_b);
}

inline real tau_stable(const real magu, const real Ri_b, const real h,
                       const real z0, const real kappa) {
  const real log_hz0 = log(h / z0);
  return (magu * magu) / (log_hz0 * log_hz0)
    * (f_tau_stable(Ri_b) / f_tau_stable((real) 0.0)) * (kappa * kappa);
}

inline real heat_flux_stable(const real ti, const real ts, const real Ri_b,
                             const real h, const real z0h, const real utau,
                             const real kappa, const real Pr) {
  return (ti - ts) / log(h / z0h)
    * (f_theta_stable(Ri_b) / fabs(f_theta_stable((real) 0.0)))
    * kappa * (utau / Pr);
}

/*
 * Similarity laws and corrections for the UNSTABLE (convective) regime:
 * Based on Louis 1979
 */

inline real f_tau_convective(const real Ri_b, const real c) {
  return (real) 1.0 - ((real) 2.0 * Ri_b) / ((real) 1.0 + c * sqrt(fabs(Ri_b)));
}

inline real f_theta_convective(const real Ri_b, const real c) {
  /* Functionally identical to f_tau_convective in Louis 1979 */
  return (real) 1.0 - ((real) 2.0 * Ri_b) / ((real) 1.0 + c * sqrt(fabs(Ri_b)));
}

inline real tau_convective(const real magu, const real Ri_b, const real h,
                           const real z0, const real kappa) {
  const real a = kappa / log(h / z0);
  const real b = (real) 2.0;
  const real c = (real) 7.4 * (a * a) * b * sqrt(h / z0);

  return (a * a) * (magu * magu) * f_tau_convective(Ri_b, c);
}

inline real heat_flux_convective(const real ti, const real ts, const real Ri_b,
                                 const real h, const real magu, const real z0h,
                                 const real kappa) {
  const real a = kappa / log(h / z0h);
  const real b = (real) 2.0;
  const real c = (real) 5.3 * (a * a) * b * sqrt(h / z0h);

  return -(a * a) / (real) 0.74 * magu * (ti - ts)
    * f_theta_convective(Ri_b, c);
}

/*
 * Similarity laws and corrections for the NEUTRAL regime:
 */

inline real tau_neutral(const real magu, const real h, const real z0,
                        const real kappa) {
  const real val = (kappa * magu) / log(h / z0);
  return val * val;
}

inline real heat_flux_neutral(const real ti, const real ts, const real h,
                              const real z0h, const real utau,
                              const real kappa) {
  return kappa * utau * (ti - ts) / log(h / z0h);
}

/**
 * Body of the Richardson wall model, shared by the two boundary
 * condition variants. @a bc_type is 0 for Neumann, 1 for Dirichlet.
 */
inline void richardson_body(__global const real * __restrict__ u_d,
                            __global const real * __restrict__ v_d,
                            __global const real * __restrict__ w_d,
                            __global const real * __restrict__ temp_d,
                            __global const real * __restrict__ temp_w_d,
                            __global const real * __restrict__ h_d,
                            __global const real * __restrict__ n_x_d,
                            __global const real * __restrict__ n_y_d,
                            __global const real * __restrict__ n_z_d,
                            __global real * __restrict__ tau_x_d,
                            __global real * __restrict__ tau_y_d,
                            __global real * __restrict__ tau_z_d,
                            const int n_nodes,
                            const real kappa,
                            __global const real * __restrict__ mu_w_d,
                            __global const real * __restrict__ rho_w_d,
                            const real g1,
                            const real g2,
                            const real g3,
                            const real Pr,
                            const real z0,
                            const real z0h_in,
                            const real bc_value,
                            __global real * __restrict__ Ri_b_diagn,
                            __global real * __restrict__ L_ob_diagn,
                            __global real * __restrict__ utau_diagn,
                            __global real * __restrict__ magu_diagn,
                            __global real * __restrict__ ti_diagn,
                            __global real * __restrict__ ts_diagn,
                            __global real * __restrict__ q_diagn,
                            const int bc_type) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  const real Ri_threshold = (real) 1e-4;

  for (int i = idx; i < n_nodes; i += str) {
    real ui = u_d[i];
    real vi = v_d[i];
    real wi = w_d[i];
    const real ti = temp_d[i];
    const real hi = h_d[i];
    const real mu = mu_w_d[i];
    const real rho = rho_w_d[i];

    /* Extract the local normal vector */
    const real nx = n_x_d[i];
    const real ny = n_y_d[i];
    const real nz = n_z_d[i];

    /* Get the tangential component */
    const real normu = ui * nx + vi * ny + wi * nz;
    ui -= normu * nx;
    vi -= normu * ny;
    wi -= normu * nz;

    real magu = sqrt(ui * ui + vi * vi + wi * wi);
    magu = fmax(magu, (real) 1e-6);

    real utau = kappa * magu / log(hi / z0);

    /* Zilitinkevich 1995 correlation for thermal roughness */
    real z0h;
    if (z0h_in < (real) 0.0) {
      /* Note that this uses previous timestep's utau, hence
         lags behind by one dt. usually very negligible */
      z0h = z0 * exp(z0h_in * sqrt((utau * z0) / (mu / rho)));
    } else {
      z0h = z0h_in;
    }

    real ts = (real) 0.0;
    real q = (real) 0.0;

    /* Initialize variables based on Boundary Condition */
    if (bc_type == 0) {         /* Neumann */
      q = bc_value;
    } else {                    /* Dirichlet */
      ts = bc_value;
      q = kappa * utau * (ts - ti) / log(hi / z0h);
    }

    real Ri_b;
    const real g_dot_n = fabs(g1 * nx + g2 * ny + g3 * nz);

    /* Compute Bulk Richardson Number */
    if (bc_type == 0) {
      Ri_b = -g_dot_n * hi / ti * q / (magu * magu * magu * kappa * kappa);
    } else {
      Ri_b = g_dot_n * hi / ti * (ti - ts) / (magu * magu);
    }

    real tau_mag = (real) 0.0;

    /* Stability regime branching */
    if (Ri_b > Ri_threshold) {          /* Stable */
      tau_mag = tau_stable(magu, Ri_b, hi, z0, kappa);
      utau = sqrt(tau_mag);
      if (bc_type == 1) {
        q = heat_flux_stable(ti, ts, Ri_b, hi, z0h, utau, kappa, Pr);
      }
    } else if (Ri_b < -Ri_threshold) {  /* Convective */
      tau_mag = tau_convective(magu, Ri_b, hi, z0, kappa);
      utau = sqrt(tau_mag);
      if (bc_type == 1) {
        q = heat_flux_convective(ti, ts, Ri_b, hi, magu, z0h, kappa);
      }
    } else {                            /* Neutral */
      tau_mag = tau_neutral(magu, hi, z0, kappa);
      utau = sqrt(tau_mag);
      if (bc_type == 1) {
        q = heat_flux_neutral(ti, ts, hi, z0h, utau, kappa);
      }
    }

    /* Apply spatial distribution */
    tau_x_d[i] = -rho * tau_mag * ui / magu;
    tau_y_d[i] = -rho * tau_mag * vi / magu;
    tau_z_d[i] = -rho * tau_mag * wi / magu;

    /* Note: L_ob is a pure diagnostic and is not computed here */
    Ri_b_diagn[i] = Ri_b;
    L_ob_diagn[i] = (real) 9999.0;
    utau_diagn[i] = utau;
    magu_diagn[i] = magu;
    ti_diagn[i] = ti;
    ts_diagn[i] = temp_w_d[i];
    q_diagn[i] = q;
  }
}

/** Device kernel for the Richardson wall model, Neumann b.c. */
__kernel void richardson_compute_neumann_kernel(
    __global const real * __restrict__ u_d,
    __global const real * __restrict__ v_d,
    __global const real * __restrict__ w_d,
    __global const real * __restrict__ temp_d,
    __global const real * __restrict__ temp_w_d,
    __global const real * __restrict__ h_d,
    __global const real * __restrict__ n_x_d,
    __global const real * __restrict__ n_y_d,
    __global const real * __restrict__ n_z_d,
    __global real * __restrict__ tau_x_d,
    __global real * __restrict__ tau_y_d,
    __global real * __restrict__ tau_z_d,
    const int n_nodes,
    const real kappa,
    __global const real * __restrict__ mu_w_d,
    __global const real * __restrict__ rho_w_d,
    const real g1, const real g2, const real g3,
    const real Pr, const real z0, const real z0h_in, const real bc_value,
    __global real * __restrict__ Ri_b_diagn,
    __global real * __restrict__ L_ob_diagn,
    __global real * __restrict__ utau_diagn,
    __global real * __restrict__ magu_diagn,
    __global real * __restrict__ ti_diagn,
    __global real * __restrict__ ts_diagn,
    __global real * __restrict__ q_diagn) {

  richardson_body(u_d, v_d, w_d, temp_d, temp_w_d, h_d, n_x_d, n_y_d, n_z_d,
                  tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, mu_w_d, rho_w_d,
                  g1, g2, g3, Pr, z0, z0h_in, bc_value,
                  Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn,
                  ti_diagn, ts_diagn, q_diagn, 0);
}

/** Device kernel for the Richardson wall model, Dirichlet b.c. */
__kernel void richardson_compute_dirichlet_kernel(
    __global const real * __restrict__ u_d,
    __global const real * __restrict__ v_d,
    __global const real * __restrict__ w_d,
    __global const real * __restrict__ temp_d,
    __global const real * __restrict__ temp_w_d,
    __global const real * __restrict__ h_d,
    __global const real * __restrict__ n_x_d,
    __global const real * __restrict__ n_y_d,
    __global const real * __restrict__ n_z_d,
    __global real * __restrict__ tau_x_d,
    __global real * __restrict__ tau_y_d,
    __global real * __restrict__ tau_z_d,
    const int n_nodes,
    const real kappa,
    __global const real * __restrict__ mu_w_d,
    __global const real * __restrict__ rho_w_d,
    const real g1, const real g2, const real g3,
    const real Pr, const real z0, const real z0h_in, const real bc_value,
    __global real * __restrict__ Ri_b_diagn,
    __global real * __restrict__ L_ob_diagn,
    __global real * __restrict__ utau_diagn,
    __global real * __restrict__ magu_diagn,
    __global real * __restrict__ ti_diagn,
    __global real * __restrict__ ts_diagn,
    __global real * __restrict__ q_diagn) {

  richardson_body(u_d, v_d, w_d, temp_d, temp_w_d, h_d, n_x_d, n_y_d, n_z_d,
                  tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, mu_w_d, rho_w_d,
                  g1, g2, g3, Pr, z0, z0h_in, bc_value,
                  Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn,
                  ti_diagn, ts_diagn, q_diagn, 1);
}

#endif // __WALL_MODELS_RICHARDSON_KERNEL_CL__
