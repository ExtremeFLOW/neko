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
 * Metal compute kernels for the MOST wall model.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Similarity laws and corrections for the STABLE regime:
 * REFERENCE: Cheng, Y., and W. Brutsaert (2005), Flux-profile relationships
 * for wind speed and temperature in the stable atmospheric boundary layer,
 * Bound.-Layer Meteorol., 3, 519-538.
 */

static float corr_m_stable(const float z, const float L_ob) {
  /* Coefficients specific to Cheng & Brutsaert (2005) */
  const float a = (float) 1.0;
  const float b = (float) 2.0 / (float) 3.0;
  const float c = (float) 5.0;
  const float d = (float) 0.35;
  const float zeta = z / L_ob;

  return -a * zeta - b * (zeta - c / d) * exp(-d * zeta) - b * c / d;
}

static float corr_h_stable(const float z, const float L_ob) {
  /* Coefficients specific to Cheng & Brutsaert (2005) */
  const float a = (float) 1.0;
  const float b = (float) 2.0 / (float) 3.0;
  const float c = (float) 5.0;
  const float d = (float) 0.35;
  const float zeta = z / L_ob;

  return -b * (zeta - c / d) * exp(-d * zeta)
    - pow((float) 1.0 + (float) 2.0 / (float) 3.0 * a * zeta, (float) 1.5)
    - b * c / d + (float) 1.0;
}

static float slaw_m_stable(const float z, const float L_ob, const float z0) {
  return log(z / z0) - corr_m_stable(z, L_ob) + corr_m_stable(z0, L_ob);
}

static float slaw_h_stable(const float z, const float L_ob, const float z0h) {
  return log(z / z0h) - corr_h_stable(z, L_ob) + corr_h_stable(z0h, L_ob);
}

static float f_neumann_stable(const float Ri_b, const float z, const float z0,
                             const float z0h, const float L_ob, const float Pr) {
  return Ri_b - Pr * z / L_ob / (slaw_m_stable(z, L_ob, z0)
                                 * slaw_m_stable(z, L_ob, z0)
                                 * slaw_m_stable(z, L_ob, z0));
}

static float dfdl_neumann_stable(const float l_upper, const float l_lower,
                                const float z, const float z0, const float z0h,
                                const float fd_h, const float Pr) {
  const float up = -z / l_upper / (slaw_m_stable(z, l_upper, z0)
                                  * slaw_m_stable(z, l_upper, z0)
                                  * slaw_m_stable(z, l_upper, z0));

  const float low = z / l_lower / (slaw_m_stable(z, l_lower, z0)
                                  * slaw_m_stable(z, l_lower, z0)
                                  * slaw_m_stable(z, l_lower, z0));

  return Pr * (up + low) / ((float) 2.0 * fd_h);
}

static float f_dirichlet_stable(const float Ri_b, const float z, const float z0,
                               const float z0h, const float L_ob,
                               const float Pr) {
  return Ri_b - Pr * z / L_ob * slaw_h_stable(z, L_ob, z0h)
    / (slaw_m_stable(z, L_ob, z0) * slaw_m_stable(z, L_ob, z0));
}

static float dfdl_dirichlet_stable(const float l_upper, const float l_lower,
                                  const float z, const float z0, const float z0h,
                                  const float fd_h, const float Pr) {
  const float up = -z / l_upper * slaw_h_stable(z, l_upper, z0h)
    / (slaw_m_stable(z, l_upper, z0) * slaw_m_stable(z, l_upper, z0));

  const float low = z / l_lower * slaw_h_stable(z, l_lower, z0h)
    / (slaw_m_stable(z, l_lower, z0) * slaw_m_stable(z, l_lower, z0));

  return Pr * (up + low) / ((float) 2.0 * fd_h);
}

/*
 * Similarity laws and corrections for the UNSTABLE (convective) regime:
 * REFERENCE: Dyer, A. J. (1974), A review of flux-profile relationships,
 * Bound.-Layer Meteorol., 7, 363-372.
 * INTEGRATION: Paulson, C. A. (1970), The mathematical representation
 * of wind speed and temperature profiles in the unstable atmospheric
 * surface layer, J. Appl. Meteorol., 9, 857-861.
 */

static float corr_m_convective(const float z, const float L_ob) {
  const float zeta = z / L_ob;
  const float pi = (float) 4.0 * atan((float) 1.0);
  /* Standard Dyer-Businger coefficient gamma = 16.0 */
  const float xi = sqrt(sqrt(((float) 1.0 - (float) 16.0 * zeta)));

  return (float) 2.0 * log((float) 0.5 * ((float) 1.0 + xi))
    + log((float) 0.5 * ((float) 1.0 + xi * xi))
    - (float) 2.0 * atan(xi) + pi / (float) 2.0;
}

static float corr_h_convective(const float z, const float L_ob) {
  const float zeta = z / L_ob;
  /* Standard Dyer-Businger coefficient gamma = 16.0 */
  const float xi = sqrt(sqrt(((float) 1.0 - (float) 16.0 * zeta)));
  return (float) 2.0 * log((float) 0.5 * ((float) 1.0 + xi * xi));
}

static float slaw_m_convective(const float z, const float L_ob, const float z0) {
  return log(z / z0) - corr_m_convective(z, L_ob)
    + corr_m_convective(z0, L_ob);
}

static float slaw_h_convective(const float z, const float L_ob, const float z0h) {
  return log(z / z0h) - corr_h_convective(z, L_ob)
    + corr_h_convective(z0h, L_ob);
}

static float f_neumann_convective(const float Ri_b, const float z, const float z0,
                                 const float z0h, const float L_ob,
                                 const float Pr) {
  return Ri_b - Pr * z / L_ob / (slaw_m_convective(z, L_ob, z0)
                                 * slaw_m_convective(z, L_ob, z0)
                                 * slaw_m_convective(z, L_ob, z0));
}

static float dfdl_neumann_convective(const float l_upper, const float l_lower,
                                    const float z, const float z0,
                                    const float z0h, const float fd_h,
                                    const float Pr) {
  const float up = -z / l_upper / (slaw_m_convective(z, l_upper, z0)
                                  * slaw_m_convective(z, l_upper, z0)
                                  * slaw_m_convective(z, l_upper, z0));

  const float low = z / l_lower / (slaw_m_convective(z, l_lower, z0)
                                  * slaw_m_convective(z, l_lower, z0)
                                  * slaw_m_convective(z, l_lower, z0));

  return Pr * (up + low) / ((float) 2.0 * fd_h);
}

static float f_dirichlet_convective(const float Ri_b, const float z,
                                   const float z0, const float z0h,
                                   const float L_ob, const float Pr) {
  return Ri_b - Pr * z / L_ob * slaw_h_convective(z, L_ob, z0h)
    / (slaw_m_convective(z, L_ob, z0) * slaw_m_convective(z, L_ob, z0));
}

static float dfdl_dirichlet_convective(const float l_upper, const float l_lower,
                                      const float z, const float z0,
                                      const float z0h, const float fd_h,
                                      const float Pr) {
  const float up = -z / l_upper * slaw_h_convective(z, l_upper, z0h)
    / (slaw_m_convective(z, l_upper, z0) * slaw_m_convective(z, l_upper, z0));

  const float low = z / l_lower * slaw_h_convective(z, l_lower, z0h)
    / (slaw_m_convective(z, l_lower, z0) * slaw_m_convective(z, l_lower, z0));

  return Pr * (up + low) / ((float) 2.0 * fd_h);
}

/*
 * Similarity laws and corrections for the NEUTRAL regime:
 */

static float slaw_m_neutral(const float z, const float z0) {
  return log(z / z0);
}

static float slaw_h_neutral(const float z, const float z0h) {
  return log(z / z0h);
}

/**
 * Body of the MOST wall model, shared by the two boundary condition
 * variants. @a bc_type is 0 for Neumann, 1 for Dirichlet.
 */
static void most_body(device const float *u_d,
                      device const float *v_d,
                      device const float *w_d,
                      device const float *temp_d,
                      device const float *temp_w_d,
                      device const float *h_d,
                      device const float *n_x_d,
                      device const float *n_y_d,
                      device const float *n_z_d,
                      device float *tau_x_d,
                      device float *tau_y_d,
                      device float *tau_z_d,
                      const int n_nodes,
                      const float kappa,
                      device const float *mu_w_d,
                      device const float *rho_w_d,
                      const float g1,
                      const float g2,
                      const float g3,
                      const float Pr,
                      const float z0,
                      const float z0h_in,
                      const float bc_value,
                      device float *Ri_b_diagn,
                      device float *L_ob_diagn,
                      device float *utau_diagn,
                      device float *magu_diagn,
                      device float *ti_diagn,
                      device float *ts_diagn,
                      device float *q_diagn,
                      const int bc_type,
                      const int i) {

  if (i >= n_nodes) return;

  const float Ri_threshold = (float) 1e-4;
  const float tol = (float) 1e-3;
  const float NR_step = (float) 1e-3;
  const int max_iter = 50;

  {
    float ui = u_d[i];
    float vi = v_d[i];
    float wi = w_d[i];
    const float ti = temp_d[i];
    const float hi = h_d[i];
    const float mu = mu_w_d[i];
    const float rho = rho_w_d[i];

    /* Extract the local normal vector */
    const float nx = n_x_d[i];
    const float ny = n_y_d[i];
    const float nz = n_z_d[i];

    /* Get the tangential component */
    const float normu = ui * nx + vi * ny + wi * nz;
    ui -= normu * nx;
    vi -= normu * ny;
    wi -= normu * nz;

    /* The magnitude used for MOST is that of the tangential vector */
    float magu = sqrt(ui * ui + vi * vi + wi * wi);
    magu = fmax(magu, (float) 1e-6);

    /* Initial utau estimate */
    float utau = kappa * magu / log(hi / z0);

    /* Zilitinkevich 1995 correlation for thermal roughness */
    float z0h;
    if (z0h_in < (float) 0.0) {
      z0h = z0 * exp(z0h_in * sqrt((utau * z0) / (mu / rho)));
    } else {
      z0h = z0h_in;
    }

    float ts = (float) 0.0;
    float q = (float) 0.0;

    if (bc_type == 0) {        /* Neumann */
      q = bc_value;
    } else {                   /* Dirichlet */
      ts = bc_value;
      q = kappa / Pr * utau * (ts - ti) / log(hi / z0h);
    }

    float Ri_b;
    const float g_dot_n = fabs(g1 * nx + g2 * ny + g3 * nz);
    if (bc_type == 0) {
      Ri_b = -g_dot_n * hi / ti * q * Pr / (magu * magu * magu * kappa * kappa);
    } else {
      Ri_b = g_dot_n * hi / ti * (ti - ts) / (magu * magu);
    }

    float L_ob = (float) 0.0;    /* neutral default */

    const float L_sign = (Ri_b > (float) 0.0) ? (float) 1.0 : (float) -1.0;

    /* Stability branching */
    if (fabs(Ri_b) <= Ri_threshold) {
      /* NEUTRAL case */
      utau = kappa * magu / slaw_m_neutral(hi, z0);
      if (bc_type == 1) {
        q = kappa / Pr * utau * (ts - ti) / slaw_h_neutral(hi, z0h);
      }
    } else {
      /* STABLE or CONVECTIVE (Newton-Raphson) */
      if (Ri_b > (float) 0.0) {
        L_ob = hi / fmax(Ri_b, Ri_threshold);
      } else {
        L_ob = hi / fmin(Ri_b, -Ri_threshold);
      }

      float L_old;
      for (int it = 0; it < max_iter; ++it) {
        L_old = L_ob;
        float f_val, dfdl;
        const float fd_h = NR_step * L_ob;
        const float L_upper = L_ob + fd_h;
        const float L_lower = L_ob - fd_h;

        /* Similarity law selected by stability and b.c. type */
        if (Ri_b > (float) 0.0) {   /* Stable */
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

        if (fabs(dfdl) < (float) 1e-12) break;
        L_ob -= f_val / dfdl;
        if (L_ob * L_sign <= (float) 0.0) L_ob = (float) 0.5 * L_old;
        L_ob = L_sign * fmax(fmin(fabs(L_ob), (float) 1e8), (float) 1e-8);
        if (fabs((L_ob - L_old) / L_ob) < tol) break;
      }

      /* Final local variables update */
      if (Ri_b > (float) 0.0) {
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


/** Metal kernel for the MOST wall model, Neumann b.c. */
kernel void most_compute_neumann_kernel(device const float *u_d [[ buffer(0) ]],
                                        device const float *v_d [[ buffer(1) ]],
                                        device const float *w_d [[ buffer(2) ]],
                                        device const float *temp_d [[ buffer(3) ]],
                                        device const float *temp_w_d [[ buffer(4) ]],
                                        device const float *h_d [[ buffer(5) ]],
                                        device const float *n_x_d [[ buffer(6) ]],
                                        device const float *n_y_d [[ buffer(7) ]],
                                        device const float *n_z_d [[ buffer(8) ]],
                                        device float *tau_x_d [[ buffer(9) ]],
                                        device float *tau_y_d [[ buffer(10) ]],
                                        device float *tau_z_d [[ buffer(11) ]],
                                        constant int &n_nodes [[ buffer(12) ]],
                                        constant float &kappa [[ buffer(13) ]],
                                        device const float *mu_w_d [[ buffer(14) ]],
                                        device const float *rho_w_d [[ buffer(15) ]],
                                        constant float &g1 [[ buffer(16) ]],
                                        constant float &g2 [[ buffer(17) ]],
                                        constant float &g3 [[ buffer(18) ]],
                                        constant float &Pr [[ buffer(19) ]],
                                        constant float &z0 [[ buffer(20) ]],
                                        constant float &z0h_in [[ buffer(21) ]],
                                        constant float &bc_value [[ buffer(22) ]],
                                        device float *Ri_b_diagn [[ buffer(23) ]],
                                        device float *L_ob_diagn [[ buffer(24) ]],
                                        device float *utau_diagn [[ buffer(25) ]],
                                        device float *magu_diagn [[ buffer(26) ]],
                                        device float *ti_diagn [[ buffer(27) ]],
                                        device float *ts_diagn [[ buffer(28) ]],
                                        device float *q_diagn [[ buffer(29) ]],
                                        uint idx [[ thread_position_in_grid ]]) {
  most_body(u_d, v_d, w_d, temp_d, temp_w_d, h_d, n_x_d, n_y_d, n_z_d,
            tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, mu_w_d, rho_w_d,
            g1, g2, g3, Pr, z0, z0h_in, bc_value,
            Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn,
            ti_diagn, ts_diagn, q_diagn, 0, (int) idx);
}

/** Metal kernel for the MOST wall model, Dirichlet b.c. */
kernel void most_compute_dirichlet_kernel(device const float *u_d [[ buffer(0) ]],
                                          device const float *v_d [[ buffer(1) ]],
                                          device const float *w_d [[ buffer(2) ]],
                                          device const float *temp_d [[ buffer(3) ]],
                                          device const float *temp_w_d [[ buffer(4) ]],
                                          device const float *h_d [[ buffer(5) ]],
                                          device const float *n_x_d [[ buffer(6) ]],
                                          device const float *n_y_d [[ buffer(7) ]],
                                          device const float *n_z_d [[ buffer(8) ]],
                                          device float *tau_x_d [[ buffer(9) ]],
                                          device float *tau_y_d [[ buffer(10) ]],
                                          device float *tau_z_d [[ buffer(11) ]],
                                          constant int &n_nodes [[ buffer(12) ]],
                                          constant float &kappa [[ buffer(13) ]],
                                          device const float *mu_w_d [[ buffer(14) ]],
                                          device const float *rho_w_d [[ buffer(15) ]],
                                          constant float &g1 [[ buffer(16) ]],
                                          constant float &g2 [[ buffer(17) ]],
                                          constant float &g3 [[ buffer(18) ]],
                                          constant float &Pr [[ buffer(19) ]],
                                          constant float &z0 [[ buffer(20) ]],
                                          constant float &z0h_in [[ buffer(21) ]],
                                          constant float &bc_value [[ buffer(22) ]],
                                          device float *Ri_b_diagn [[ buffer(23) ]],
                                          device float *L_ob_diagn [[ buffer(24) ]],
                                          device float *utau_diagn [[ buffer(25) ]],
                                          device float *magu_diagn [[ buffer(26) ]],
                                          device float *ti_diagn [[ buffer(27) ]],
                                          device float *ts_diagn [[ buffer(28) ]],
                                          device float *q_diagn [[ buffer(29) ]],
                                          uint idx [[ thread_position_in_grid ]]) {
  most_body(u_d, v_d, w_d, temp_d, temp_w_d, h_d, n_x_d, n_y_d, n_z_d,
            tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, mu_w_d, rho_w_d,
            g1, g2, g3, Pr, z0, z0h_in, bc_value,
            Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn,
            ti_diagn, ts_diagn, q_diagn, 1, (int) idx);
}
