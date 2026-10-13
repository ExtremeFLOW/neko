#ifndef __LES_VREMAN_NUT_KERNEL_CL__
#define __LES_VREMAN_NUT_KERNEL_CL__
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
 * Device kernel for vreman_nut_compute
 */
__kernel void vreman_nut_compute_kernel(__global const real * __restrict__ a11,
                                        __global const real * __restrict__ a12,
                                        __global const real * __restrict__ a13,
                                        __global const real * __restrict__ a21,
                                        __global const real * __restrict__ a22,
                                        __global const real * __restrict__ a23,
                                        __global const real * __restrict__ a31,
                                        __global const real * __restrict__ a32,
                                        __global const real * __restrict__ a33,
                                        __global const real * __restrict__ delta,
                                        __global real * __restrict__ nut,
                                        __global const real * __restrict__ mult,
                                        const real c,
                                        const real eps,
                                        const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  const real zero = (real) 0.0;

  for (int i = idx; i < n; i += str) {

    const real a11_r = a11[i];
    const real a12_r = a12[i];
    const real a13_r = a13[i];
    const real a21_r = a21[i];
    const real a22_r = a22[i];
    const real a23_r = a23[i];
    const real a31_r = a31[i];
    const real a32_r = a32[i];
    const real a33_r = a33[i];
    const real delta_r = delta[i];

    const real beta11 = a11_r * a11_r + a21_r * a21_r + a31_r * a31_r;
    const real beta22 = a12_r * a12_r + a22_r * a22_r + a32_r * a32_r;
    const real beta33 = a13_r * a13_r + a23_r * a23_r + a33_r * a33_r;
    const real beta12 = a11_r * a12_r + a21_r * a22_r + a31_r * a32_r;
    const real beta13 = a11_r * a13_r + a21_r * a23_r + a31_r * a33_r;
    const real beta23 = a12_r * a13_r + a22_r * a23_r + a32_r * a33_r;

    real b_beta = beta11 * beta22 - beta12 * beta12
      + beta11 * beta33 - beta13 * beta13
      + beta22 * beta33 - beta23 * beta23;

    b_beta = max(zero, b_beta);

    const real aijaij = beta11 + beta22 + beta33;

    nut[i] = c * delta_r * delta_r
      * sqrt(b_beta / (aijaij + eps)) * mult[i];
  }
}

/**
 * Device kernel for vreman_nut_compute_buoy
 */
__kernel void vreman_nut_compute_buoy_kernel(__global const real * __restrict__ a11,
                                             __global const real * __restrict__ a12,
                                             __global const real * __restrict__ a13,
                                             __global const real * __restrict__ a21,
                                             __global const real * __restrict__ a22,
                                             __global const real * __restrict__ a23,
                                             __global const real * __restrict__ a31,
                                             __global const real * __restrict__ a32,
                                             __global const real * __restrict__ a33,
                                             __global const real * __restrict__ delta,
                                             __global real * __restrict__ nut,
                                             __global const real * __restrict__ mult,
                                             const real c,
                                             const real eps,
                                             const int n,
                                             __global const real * __restrict__ dTdx,
                                             __global const real * __restrict__ dTdy,
                                             __global const real * __restrict__ dTdz,
                                             const real n1,
                                             const real n2,
                                             const real n3,
                                             const real g1,
                                             const real g2,
                                             const real g3,
                                             const real ri_c,
                                             const real ref_temp) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  const real zero = (real) 0.0;
  const real one = (real) 1.0;

  for (int i = idx; i < n; i += str) {

    const real a11_r = a11[i];
    const real a12_r = a12[i];
    const real a13_r = a13[i];
    const real a21_r = a21[i];
    const real a22_r = a22[i];
    const real a23_r = a23[i];
    const real a31_r = a31[i];
    const real a32_r = a32[i];
    const real a33_r = a33[i];
    const real delta_r = delta[i];

    const real beta11 = a11_r * a11_r + a21_r * a21_r + a31_r * a31_r;
    const real beta22 = a12_r * a12_r + a22_r * a22_r + a32_r * a32_r;
    const real beta33 = a13_r * a13_r + a23_r * a23_r + a33_r * a33_r;
    const real beta12 = a11_r * a12_r + a21_r * a22_r + a31_r * a32_r;
    const real beta13 = a11_r * a13_r + a21_r * a23_r + a31_r * a33_r;
    const real beta23 = a12_r * a13_r + a22_r * a23_r + a32_r * a33_r;

    real b_beta = beta11 * beta22 - beta12 * beta12
      + beta11 * beta33 - beta13 * beta13
      + beta22 * beta33 - beta23 * beta23;

    b_beta = max(zero, b_beta);

    const real aijaij = beta11 + beta22 + beta33;

    const real nut0 = c * delta_r * delta_r
      * sqrt(b_beta / (aijaij + eps)) * mult[i];

    /* Scalar gradient for buoyancy */
    const real buoyancy = (g1 * dTdx[i] + g2 * dTdy[i] + g3 * dTdz[i])
      / ref_temp;

    /* Directional derivative along n */
    const real du_n1 = a11_r * n1 + a12_r * n2 + a13_r * n3;
    const real du_n2 = a21_r * n1 + a22_r * n2 + a23_r * n3;
    const real du_n3 = a31_r * n1 + a32_r * n2 + a33_r * n3;

    const real du_parallel = du_n1 * n1 + du_n2 * n2 + du_n3 * n3;

    const real sh1 = du_n1 - du_parallel * n1;
    const real sh2 = du_n2 - du_parallel * n2;
    const real sh3 = du_n3 - du_parallel * n3;

    const real shear_sq = sh1 * sh1 + sh2 * sh2 + sh3 * sh3;

    const real ri = buoyancy / (shear_sq + eps);

    real out;
    if (ri <= ri_c) {
      real v = one - ri / ri_c;
      if (v < zero) v = zero;
      out = nut0 * sqrt(v);
    } else {
      out = eps;
    }

    nut[i] = out;
  }
}

#endif // __LES_VREMAN_NUT_KERNEL_CL__
