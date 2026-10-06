#ifndef __LES_DYNAMIC_SMAGORINSKY_NUT_KERNEL_CL__
#define __LES_DYNAMIC_SMAGORINSKY_NUT_KERNEL_CL__
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
 * Device kernel for s_abs_compute
 */
__kernel void s_abs_compute_kernel(__global real * __restrict__ s_abs,
                                   __global const real * __restrict__ s11,
                                   __global const real * __restrict__ s22,
                                   __global const real * __restrict__ s33,
                                   __global const real * __restrict__ s12,
                                   __global const real * __restrict__ s13,
                                   __global const real * __restrict__ s23,
                                   const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    const real s11_r = s11[i];
    const real s22_r = s22[i];
    const real s33_r = s33[i];
    const real s12_r = s12[i];
    const real s13_r = s13[i];
    const real s23_r = s23[i];

    s_abs[i] = sqrt((real) 2.0 * (s11_r * s11_r +
                                  s22_r * s22_r +
                                  s33_r * s33_r) +
                    (real) 4.0 * (s12_r * s12_r +
                                  s13_r * s13_r +
                                  s23_r * s23_r));
  }
}

/**
 * Device kernel for lij_compute_part1
 */
__kernel void lij_compute_part1_kernel(__global real * __restrict__ l11,
                                       __global real * __restrict__ l22,
                                       __global real * __restrict__ l33,
                                       __global real * __restrict__ l12,
                                       __global real * __restrict__ l13,
                                       __global real * __restrict__ l23,
                                       __global const real * __restrict__ u,
                                       __global const real * __restrict__ v,
                                       __global const real * __restrict__ w,
                                       __global const real * __restrict__ fu,
                                       __global const real * __restrict__ fv,
                                       __global const real * __restrict__ fw,
                                       __global real * __restrict__ fuu,
                                       __global real * __restrict__ fvv,
                                       __global real * __restrict__ fww,
                                       __global real * __restrict__ fuv,
                                       __global real * __restrict__ fuw,
                                       __global real * __restrict__ fvw,
                                       const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    const real u_r = u[i];
    const real v_r = v[i];
    const real w_r = w[i];
    const real fu_r = fu[i];
    const real fv_r = fv[i];
    const real fw_r = fw[i];

    l11[i] = fu_r * fu_r;
    l22[i] = fv_r * fv_r;
    l33[i] = fw_r * fw_r;
    l12[i] = fu_r * fv_r;
    l13[i] = fu_r * fw_r;
    l23[i] = fv_r * fw_r;

    fuu[i] = u_r * u_r;
    fvv[i] = v_r * v_r;
    fww[i] = w_r * w_r;
    fuv[i] = u_r * v_r;
    fuw[i] = u_r * w_r;
    fvw[i] = v_r * w_r;
  }
}

/**
 * Device kernel for lij_compute_part2
 */
__kernel void lij_compute_part2_kernel(__global real * __restrict__ l11,
                                       __global real * __restrict__ l22,
                                       __global real * __restrict__ l33,
                                       __global real * __restrict__ l12,
                                       __global real * __restrict__ l13,
                                       __global real * __restrict__ l23,
                                       __global const real * __restrict__ fuu,
                                       __global const real * __restrict__ fvv,
                                       __global const real * __restrict__ fww,
                                       __global const real * __restrict__ fuv,
                                       __global const real * __restrict__ fuw,
                                       __global const real * __restrict__ fvw,
                                       const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    l11[i] -= fuu[i];
    l22[i] -= fvv[i];
    l33[i] -= fww[i];
    l12[i] -= fuv[i];
    l13[i] -= fuw[i];
    l23[i] -= fvw[i];
  }
}

/**
 * Device kernel for mij_compute_part1
 */
__kernel void mij_compute_part1_kernel(__global real * __restrict__ m11,
                                       __global real * __restrict__ m22,
                                       __global real * __restrict__ m33,
                                       __global real * __restrict__ m12,
                                       __global real * __restrict__ m13,
                                       __global real * __restrict__ m23,
                                       __global const real * __restrict__ s_abs,
                                       __global const real * __restrict__ s11,
                                       __global const real * __restrict__ s22,
                                       __global const real * __restrict__ s33,
                                       __global const real * __restrict__ s12,
                                       __global const real * __restrict__ s13,
                                       __global const real * __restrict__ s23,
                                       __global const real * __restrict__ fs_abs,
                                       __global const real * __restrict__ fs11,
                                       __global const real * __restrict__ fs22,
                                       __global const real * __restrict__ fs33,
                                       __global const real * __restrict__ fs12,
                                       __global const real * __restrict__ fs13,
                                       __global const real * __restrict__ fs23,
                                       __global real * __restrict__ fsabss11,
                                       __global real * __restrict__ fsabss22,
                                       __global real * __restrict__ fsabss33,
                                       __global real * __restrict__ fsabss12,
                                       __global real * __restrict__ fsabss13,
                                       __global real * __restrict__ fsabss23,
                                       const real delta_ratio2,
                                       const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    const real s_abs_r = s_abs[i];
    const real fs_abs_r = fs_abs[i];

    m11[i] = delta_ratio2 * (fs_abs_r * fs11[i]);
    m22[i] = delta_ratio2 * (fs_abs_r * fs22[i]);
    m33[i] = delta_ratio2 * (fs_abs_r * fs33[i]);
    m12[i] = delta_ratio2 * (fs_abs_r * fs12[i]);
    m13[i] = delta_ratio2 * (fs_abs_r * fs13[i]);
    m23[i] = delta_ratio2 * (fs_abs_r * fs23[i]);

    fsabss11[i] = s_abs_r * s11[i];
    fsabss22[i] = s_abs_r * s22[i];
    fsabss33[i] = s_abs_r * s33[i];
    fsabss12[i] = s_abs_r * s12[i];
    fsabss13[i] = s_abs_r * s13[i];
    fsabss23[i] = s_abs_r * s23[i];
  }
}

/**
 * Device kernel for mij_nut_compute_part2
 */
__kernel void mij_nut_compute_part2_kernel(__global const real * __restrict__ m11,
                                           __global const real * __restrict__ m22,
                                           __global const real * __restrict__ m33,
                                           __global const real * __restrict__ m12,
                                           __global const real * __restrict__ m13,
                                           __global const real * __restrict__ m23,
                                           __global const real * __restrict__ l11,
                                           __global const real * __restrict__ l22,
                                           __global const real * __restrict__ l33,
                                           __global const real * __restrict__ l12,
                                           __global const real * __restrict__ l13,
                                           __global const real * __restrict__ l23,
                                           __global const real * __restrict__ fsabss11,
                                           __global const real * __restrict__ fsabss22,
                                           __global const real * __restrict__ fsabss33,
                                           __global const real * __restrict__ fsabss12,
                                           __global const real * __restrict__ fsabss13,
                                           __global const real * __restrict__ fsabss23,
                                           __global real * __restrict__ num,
                                           __global real * __restrict__ den,
                                           __global real * __restrict__ c_dyn,
                                           __global const real * __restrict__ delta,
                                           __global const real * __restrict__ s_abs,
                                           __global real * __restrict__ nut,
                                           const real alpha,
                                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  const real zero = (real) 0.0;
  const real one = (real) 1.0;
  const real two = (real) 2.0;

  for (int i = idx; i < n; i += str) {
    const real delta_r = delta[i];
    const real delta2 = delta_r * delta_r;
    real c_dyn_r;

    const real m11_r = delta2 * (m11[i] - fsabss11[i]);
    const real m22_r = delta2 * (m22[i] - fsabss22[i]);
    const real m33_r = delta2 * (m33[i] - fsabss33[i]);
    const real m12_r = delta2 * (m12[i] - fsabss12[i]);
    const real m13_r = delta2 * (m13[i] - fsabss13[i]);
    const real m23_r = delta2 * (m23[i] - fsabss23[i]);

    const real num_curr = m11_r * l11[i]
      + m22_r * l22[i]
      + m33_r * l33[i]
      + two * (m12_r * l12[i]
               + m13_r * l13[i]
               + m23_r * l23[i]);

    const real den_curr = m11_r * m11_r
      + m22_r * m22_r
      + m33_r * m33_r
      + two * (m12_r * m12_r
               + m13_r * m13_r
               + m23_r * m23_r);

    const real num_r = alpha * num[i] + (one - alpha) * num_curr;
    const real den_r = alpha * den[i] + (one - alpha) * den_curr;

    num[i] = num_r;
    den[i] = den_r;

    if (den_r > zero) {
      c_dyn_r = (real) 0.5 * num_r / den_r;
    } else {
      c_dyn_r = zero;
    }

    c_dyn_r = max(c_dyn_r, zero);
    c_dyn[i] = c_dyn_r;
    nut[i] = c_dyn_r * delta2 * s_abs[i];
  }
}

#endif // __LES_DYNAMIC_SMAGORINSKY_NUT_KERNEL_CL__
