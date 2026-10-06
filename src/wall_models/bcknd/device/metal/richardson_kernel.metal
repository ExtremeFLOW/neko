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
 * Metal compute kernels for the Richardson wall model.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Similarity laws and corrections for the STABLE regime:
 * Based on Mauritsen et al. 2007
 */

static float f_tau_stable(const float Ri_b) {
  return (float) 0.17 * ((float) 0.25 +
                        (float) 0.75 / ((float) 1.0 + (float) 4.0 * Ri_b));
}

static float f_theta_stable(const float Ri_b) {
  return (float) -0.145 / ((float) 1.0 + (float) 4.0 * Ri_b);
}

static float tau_stable(const float magu, const float Ri_b, const float h,
                       const float z0, const float kappa) {
  const float log_hz0 = log(h / z0);
  return (magu * magu) / (log_hz0 * log_hz0)
    * (f_tau_stable(Ri_b) / f_tau_stable((float) 0.0)) * (kappa * kappa);
}

static float heat_flux_stable(const float ti, const float ts, const float Ri_b,
                             const float h, const float z0h, const float utau,
                             const float kappa, const float Pr) {
  return (ti - ts) / log(h / z0h)
    * (f_theta_stable(Ri_b) / fabs(f_theta_stable((float) 0.0)))
    * kappa * (utau / Pr);
}

/*
 * Similarity laws and corrections for the UNSTABLE (convective) regime:
 * Based on Louis 1979
 */

static float f_tau_convective(const float Ri_b, const float c) {
  return (float) 1.0 - ((float) 2.0 * Ri_b) / ((float) 1.0 + c * sqrt(fabs(Ri_b)));
}

static float f_theta_convective(const float Ri_b, const float c) {
  /* Functionally identical to f_tau_convective in Louis 1979 */
  return (float) 1.0 - ((float) 2.0 * Ri_b) / ((float) 1.0 + c * sqrt(fabs(Ri_b)));
}

static float tau_convective(const float magu, const float Ri_b, const float h,
                           const float z0, const float kappa) {
  const float a = kappa / log(h / z0);
  const float b = (float) 2.0;
  const float c = (float) 7.4 * (a * a) * b * sqrt(h / z0);

  return (a * a) * (magu * magu) * f_tau_convective(Ri_b, c);
}

static float heat_flux_convective(const float ti, const float ts, const float Ri_b,
                                 const float h, const float magu, const float z0h,
                                 const float kappa) {
  const float a = kappa / log(h / z0h);
  const float b = (float) 2.0;
  const float c = (float) 5.3 * (a * a) * b * sqrt(h / z0h);

  return -(a * a) / (float) 0.74 * magu * (ti - ts)
    * f_theta_convective(Ri_b, c);
}

/*
 * Similarity laws and corrections for the NEUTRAL regime:
 */

static float tau_neutral(const float magu, const float h, const float z0,
                        const float kappa) {
  const float val = (kappa * magu) / log(h / z0);
  return val * val;
}

static float heat_flux_neutral(const float ti, const float ts, const float h,
                              const float z0h, const float utau,
                              const float kappa) {
  return kappa * utau * (ti - ts) / log(h / z0h);
}

/**
 * Body of the Richardson wall model, shared by the two boundary
 * condition variants. @a bc_type is 0 for Neumann, 1 for Dirichlet.
 */
static void richardson_body(device const float *u_d,
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

    float magu = sqrt(ui * ui + vi * vi + wi * wi);
    magu = fmax(magu, (float) 1e-6);

    float utau = kappa * magu / log(hi / z0);

    /* Zilitinkevich 1995 correlation for thermal roughness */
    float z0h;
    if (z0h_in < (float) 0.0) {
      /* Note that this uses previous timestep's utau, hence
         lags behind by one dt. usually very negligible */
      z0h = z0 * exp(z0h_in * sqrt((utau * z0) / (mu / rho)));
    } else {
      z0h = z0h_in;
    }

    float ts = (float) 0.0;
    float q = (float) 0.0;

    /* Initialize variables based on Boundary Condition */
    if (bc_type == 0) {         /* Neumann */
      q = bc_value;
    } else {                    /* Dirichlet */
      ts = bc_value;
      q = kappa * utau * (ts - ti) / log(hi / z0h);
    }

    float Ri_b;
    const float g_dot_n = fabs(g1 * nx + g2 * ny + g3 * nz);

    /* Compute Bulk Richardson Number */
    if (bc_type == 0) {
      Ri_b = -g_dot_n * hi / ti * q / (magu * magu * magu * kappa * kappa);
    } else {
      Ri_b = g_dot_n * hi / ti * (ti - ts) / (magu * magu);
    }

    float tau_mag = (float) 0.0;

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
    L_ob_diagn[i] = (float) 9999.0;
    utau_diagn[i] = utau;
    magu_diagn[i] = magu;
    ti_diagn[i] = ti;
    ts_diagn[i] = temp_w_d[i];
    q_diagn[i] = q;
  }
}


/** Metal kernel for the Richardson wall model, Neumann b.c. */
kernel void richardson_compute_neumann_kernel(device const float *u_d [[ buffer(0) ]],
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
  richardson_body(u_d, v_d, w_d, temp_d, temp_w_d, h_d, n_x_d, n_y_d, n_z_d,
            tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, mu_w_d, rho_w_d,
            g1, g2, g3, Pr, z0, z0h_in, bc_value,
            Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn,
            ti_diagn, ts_diagn, q_diagn, 0, (int) idx);
}

/** Metal kernel for the Richardson wall model, Dirichlet b.c. */
kernel void richardson_compute_dirichlet_kernel(device const float *u_d [[ buffer(0) ]],
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
  richardson_body(u_d, v_d, w_d, temp_d, temp_w_d, h_d, n_x_d, n_y_d, n_z_d,
            tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, mu_w_d, rho_w_d,
            g1, g2, g3, Pr, z0, z0h_in, bc_value,
            Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn,
            ti_diagn, ts_diagn, q_diagn, 1, (int) idx);
}
