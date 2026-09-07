#ifndef __WALL_MODELS_MOST_KERNEL_CL__
#define __WALL_MODELS_MOST_KERNEL_CL__
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
 * REFERENCE: Cheng, Y., and W. Brutsaert (2005), Flux-profile relationships
 * for wind speed and temperature in the stable atmospheric boundary layer,
 * Bound.-Layer Meteorol., 3, 519-538.
 */

inline real corr_m_stable(const real z, const real L_ob) {
  /* Coefficients specific to Cheng & Brutsaert (2005) */
  const real a = (real) 1.0;
  const real b = (real) 2.0 / (real) 3.0;
  const real c = (real) 5.0;
  const real d = (real) 0.35;
  const real zeta = z / L_ob;

  return -a * zeta - b * (zeta - c / d) * exp(-d * zeta) - b * c / d;
}

inline real corr_h_stable(const real z, const real L_ob) {
  /* Coefficients specific to Cheng & Brutsaert (2005) */
  const real a = (real) 1.0;
  const real b = (real) 2.0 / (real) 3.0;
  const real c = (real) 5.0;
  const real d = (real) 0.35;
  const real zeta = z / L_ob;

  return -b * (zeta - c / d) * exp(-d * zeta)
    - pow((real) 1.0 + (real) 2.0 / (real) 3.0 * a * zeta, (real) 1.5)
    - b * c / d + (real) 1.0;
}

inline real slaw_m_stable(const real z, const real L_ob, const real z0) {
  return log(z / z0) - corr_m_stable(z, L_ob) + corr_m_stable(z0, L_ob);
}

inline real slaw_h_stable(const real z, const real L_ob, const real z0h) {
  return log(z / z0h) - corr_h_stable(z, L_ob) + corr_h_stable(z0h, L_ob);
}

inline real f_neumann_stable(const real Ri_b, const real z, const real z0,
                             const real z0h, const real L_ob, const real Pr) {
  return Ri_b - Pr * z / L_ob / (slaw_m_stable(z, L_ob, z0)
                                 * slaw_m_stable(z, L_ob, z0)
                                 * slaw_m_stable(z, L_ob, z0));
}

inline real dfdl_neumann_stable(const real l_upper, const real l_lower,
                                const real z, const real z0, const real z0h,
                                const real fd_h, const real Pr) {
  const real up = -z / l_upper / (slaw_m_stable(z, l_upper, z0)
                                  * slaw_m_stable(z, l_upper, z0)
                                  * slaw_m_stable(z, l_upper, z0));

  const real low = z / l_lower / (slaw_m_stable(z, l_lower, z0)
                                  * slaw_m_stable(z, l_lower, z0)
                                  * slaw_m_stable(z, l_lower, z0));

  return Pr * (up + low) / ((real) 2.0 * fd_h);
}

inline real f_dirichlet_stable(const real Ri_b, const real z, const real z0,
                               const real z0h, const real L_ob,
                               const real Pr) {
  return Ri_b - Pr * z / L_ob * slaw_h_stable(z, L_ob, z0h)
    / (slaw_m_stable(z, L_ob, z0) * slaw_m_stable(z, L_ob, z0));
}

inline real dfdl_dirichlet_stable(const real l_upper, const real l_lower,
                                  const real z, const real z0, const real z0h,
                                  const real fd_h, const real Pr) {
  const real up = -z / l_upper * slaw_h_stable(z, l_upper, z0h)
    / (slaw_m_stable(z, l_upper, z0) * slaw_m_stable(z, l_upper, z0));

  const real low = z / l_lower * slaw_h_stable(z, l_lower, z0h)
    / (slaw_m_stable(z, l_lower, z0) * slaw_m_stable(z, l_lower, z0));

  return Pr * (up + low) / ((real) 2.0 * fd_h);
}

/*
 * Similarity laws and corrections for the UNSTABLE (convective) regime:
 * REFERENCE: Dyer, A. J. (1974), A review of flux-profile relationships,
 * Bound.-Layer Meteorol., 7, 363-372.
 * INTEGRATION: Paulson, C. A. (1970), The mathematical representation
 * of wind speed and temperature profiles in the unstable atmospheric
 * surface layer, J. Appl. Meteorol., 9, 857-861.
 */

inline real corr_m_convective(const real z, const real L_ob) {
  const real zeta = z / L_ob;
  const real pi = (real) 4.0 * atan((real) 1.0);
  /* Standard Dyer-Businger coefficient gamma = 16.0 */
  const real xi = sqrt(sqrt(((real) 1.0 - (real) 16.0 * zeta)));

  return (real) 2.0 * log((real) 0.5 * ((real) 1.0 + xi))
    + log((real) 0.5 * ((real) 1.0 + xi * xi))
    - (real) 2.0 * atan(xi) + pi / (real) 2.0;
}

inline real corr_h_convective(const real z, const real L_ob) {
  const real zeta = z / L_ob;
  /* Standard Dyer-Businger coefficient gamma = 16.0 */
  const real xi = sqrt(sqrt(((real) 1.0 - (real) 16.0 * zeta)));
  return (real) 2.0 * log((real) 0.5 * ((real) 1.0 + xi * xi));
}

inline real slaw_m_convective(const real z, const real L_ob, const real z0) {
  return log(z / z0) - corr_m_convective(z, L_ob)
    + corr_m_convective(z0, L_ob);
}

inline real slaw_h_convective(const real z, const real L_ob, const real z0h) {
  return log(z / z0h) - corr_h_convective(z, L_ob)
    + corr_h_convective(z0h, L_ob);
}

inline real f_neumann_convective(const real Ri_b, const real z, const real z0,
                                 const real z0h, const real L_ob,
                                 const real Pr) {
  return Ri_b - Pr * z / L_ob / (slaw_m_convective(z, L_ob, z0)
                                 * slaw_m_convective(z, L_ob, z0)
                                 * slaw_m_convective(z, L_ob, z0));
}

inline real dfdl_neumann_convective(const real l_upper, const real l_lower,
                                    const real z, const real z0,
                                    const real z0h, const real fd_h,
                                    const real Pr) {
  const real up = -z / l_upper / (slaw_m_convective(z, l_upper, z0)
                                  * slaw_m_convective(z, l_upper, z0)
                                  * slaw_m_convective(z, l_upper, z0));

  const real low = z / l_lower / (slaw_m_convective(z, l_lower, z0)
                                  * slaw_m_convective(z, l_lower, z0)
                                  * slaw_m_convective(z, l_lower, z0));

  return Pr * (up + low) / ((real) 2.0 * fd_h);
}

inline real f_dirichlet_convective(const real Ri_b, const real z,
                                   const real z0, const real z0h,
                                   const real L_ob, const real Pr) {
  return Ri_b - Pr * z / L_ob * slaw_h_convective(z, L_ob, z0h)
    / (slaw_m_convective(z, L_ob, z0) * slaw_m_convective(z, L_ob, z0));
}

inline real dfdl_dirichlet_convective(const real l_upper, const real l_lower,
                                      const real z, const real z0,
                                      const real z0h, const real fd_h,
                                      const real Pr) {
  const real up = -z / l_upper * slaw_h_convective(z, l_upper, z0h)
    / (slaw_m_convective(z, l_upper, z0) * slaw_m_convective(z, l_upper, z0));

  const real low = z / l_lower * slaw_h_convective(z, l_lower, z0h)
    / (slaw_m_convective(z, l_lower, z0) * slaw_m_convective(z, l_lower, z0));

  return Pr * (up + low) / ((real) 2.0 * fd_h);
}

/*
 * Similarity laws and corrections for the NEUTRAL regime:
 */

inline real slaw_m_neutral(const real z, const real z0) {
  return log(z / z0);
}

inline real slaw_h_neutral(const real z, const real z0h) {
  return log(z / z0h);
}

/**
 * Body of the MOST wall model, shared by the two boundary condition
 * variants. @a bc_type is 0 for Neumann, 1 for Dirichlet.
 */
inline void most_body(__global const real * __restrict__ u_d,
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
  const real tol = (real) 1e-3;
  const real NR_step = (real) 1e-3;
  const int max_iter = 50;

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

    /* The magnitude used for MOST is that of the tangential vector */
    real magu = sqrt(ui * ui + vi * vi + wi * wi);
    magu = fmax(magu, (real) 1e-6);

    /* Initial utau estimate */
    real utau = kappa * magu / log(hi / z0);

    /* Zilitinkevich 1995 correlation for thermal roughness */
    real z0h;
    if (z0h_in < (real) 0.0) {
      z0h = z0 * exp(z0h_in * sqrt((utau * z0) / (mu / rho)));
    } else {
      z0h = z0h_in;
    }

    real ts = (real) 0.0;
    real q = (real) 0.0;

    if (bc_type == 0) {        /* Neumann */
      q = bc_value;
    } else {                   /* Dirichlet */
      ts = bc_value;
      q = kappa / Pr * utau * (ts - ti) / log(hi / z0h);
    }

    real Ri_b;
    const real g_dot_n = fabs(g1 * nx + g2 * ny + g3 * nz);
    if (bc_type == 0) {
      Ri_b = -g_dot_n * hi / ti * q * Pr / (magu * magu * magu * kappa * kappa);
    } else {
      Ri_b = g_dot_n * hi / ti * (ti - ts) / (magu * magu);
    }

    real L_ob = (real) 0.0;    /* neutral default */

    const real L_sign = (Ri_b > (real) 0.0) ? (real) 1.0 : (real) -1.0;

    /* Stability branching */
    if (fabs(Ri_b) <= Ri_threshold) {
      /* NEUTRAL case */
      utau = kappa * magu / slaw_m_neutral(hi, z0);
      if (bc_type == 1) {
        q = kappa / Pr * utau * (ts - ti) / slaw_h_neutral(hi, z0h);
      }
    } else {
      /* STABLE or CONVECTIVE (Newton-Raphson) */
      if (Ri_b > (real) 0.0) {
        L_ob = hi / fmax(Ri_b, Ri_threshold);
      } else {
        L_ob = hi / fmin(Ri_b, -Ri_threshold);
      }

      real L_old;
      for (int it = 0; it < max_iter; ++it) {
        L_old = L_ob;
        real f_val, dfdl;
        const real fd_h = NR_step * L_ob;
        const real L_upper = L_ob + fd_h;
        const real L_lower = L_ob - fd_h;

        /* Similarity law selected by stability and b.c. type */
        if (Ri_b > (real) 0.0) {   /* Stable */
          if (bc_type == 0) {
            f_val = f_neumann_stable(Ri_b, hi, z0, z0h, L_ob, Pr);
            dfdl = dfdl_neumann_stable(L_upper, L_lower, hi, z0, z0h,
                                       fd_h, Pr);
          } else {
            f_val = f_dirichlet_stable(Ri_b, hi, z0, z0h, L_ob, Pr);
            dfdl = dfdl_dirichlet_stable(L_upper, L_lower, hi, z0, z0h,
                                         fd_h, Pr);
          }
        } else {                   /* Convective */
          if (bc_type == 0) {
            f_val = f_neumann_convective(Ri_b, hi, z0, z0h, L_ob, Pr);
            dfdl = dfdl_neumann_convective(L_upper, L_lower, hi, z0, z0h,
                                           fd_h, Pr);
          } else {
            f_val = f_dirichlet_convective(Ri_b, hi, z0, z0h, L_ob, Pr);
            dfdl = dfdl_dirichlet_convective(L_upper, L_lower, hi, z0, z0h,
                                             fd_h, Pr);
          }
        }

        if (fabs(dfdl) < (real) 1e-12) break;
        L_ob -= f_val / dfdl;
        if (L_ob * L_sign <= (real) 0.0) L_ob = (real) 0.5 * L_old;
        L_ob = L_sign * fmax(fmin(fabs(L_ob), (real) 1e8), (real) 1e-8);
        if (fabs((L_ob - L_old) / L_ob) < tol) break;
      }

      /* Final local variables update */
      if (Ri_b > (real) 0.0) {
        utau = kappa * magu / slaw_m_stable(hi, L_ob, z0);
        if (bc_type == 1) {
          q = kappa / Pr * utau * (ts - ti) / slaw_h_stable(hi, L_ob, z0h);
        }
      } else {
        utau = kappa * magu / slaw_m_convective(hi, L_ob, z0);
        if (bc_type == 1) {
          q = kappa / Pr * utau * (ts - ti) / slaw_h_convective(hi, L_ob, z0h);
        }
      }
    }

    tau_x_d[i] = -rho * utau * utau * ui / magu;
    tau_y_d[i] = -rho * utau * utau * vi / magu;
    tau_z_d[i] = -rho * utau * utau * wi / magu;

    Ri_b_diagn[i] = Ri_b;
    L_ob_diagn[i] = L_ob;
    utau_diagn[i] = utau;
    magu_diagn[i] = magu;
    ti_diagn[i] = ti;
    ts_diagn[i] = temp_w_d[i];
    q_diagn[i] = q;
  }
}

/** Device kernel for the MOST wall model, Neumann b.c. */
__kernel void most_compute_neumann_kernel(
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

  most_body(u_d, v_d, w_d, temp_d, temp_w_d, h_d, n_x_d, n_y_d, n_z_d,
            tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, mu_w_d, rho_w_d,
            g1, g2, g3, Pr, z0, z0h_in, bc_value,
            Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn,
            ti_diagn, ts_diagn, q_diagn, 0);
}

/** Device kernel for the MOST wall model, Dirichlet b.c. */
__kernel void most_compute_dirichlet_kernel(
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

  most_body(u_d, v_d, w_d, temp_d, temp_w_d, h_d, n_x_d, n_y_d, n_z_d,
            tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, mu_w_d, rho_w_d,
            g1, g2, g3, Pr, z0, z0h_in, bc_value,
            Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn,
            ti_diagn, ts_diagn, q_diagn, 1);
}

#endif // __WALL_MODELS_MOST_KERNEL_CL__
