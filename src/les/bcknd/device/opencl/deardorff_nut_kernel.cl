#ifndef __LES_DEARDORFF_NUT_KERNEL_CL__
#define __LES_DEARDORFF_NUT_KERNEL_CL__
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
 * Device kernel for deardorff_nut_compute
 */
__kernel void deardorff_nut_compute_kernel(__global real * __restrict__ TKE,
                                           __global const real * __restrict__ dTdx,
                                           __global const real * __restrict__ dTdy,
                                           __global const real * __restrict__ dTdz,
                                           __global const real * __restrict__ a11,
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
                                           __global real * __restrict__ temperature_alphat,
                                           __global real * __restrict__ TKE_alphat,
                                           __global real * __restrict__ TKE_source,
                                           const real c_k,
                                           const real T0,
                                           const real g1,
                                           const real g2,
                                           const real g3,
                                           const real eps,
                                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  const real zero = (real) 0.0;
  const real one = (real) 1.0;
  const real two = (real) 2.0;

  for (int i = idx; i < n; i += str) {
    const real dTdx_r = dTdx[i];
    const real dTdy_r = dTdy[i];
    const real dTdz_r = dTdz[i];
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

    if (TKE[i] < eps) {
      TKE[i] = eps;
    }
    const real TKE_r = TKE[i];

    const real N2 = (dTdx_r * g1 + dTdy_r * g2 + dTdz_r * g3) / T0;
    real l;
    if (N2 > zero) {
      l = (real) 0.76 * sqrt(TKE_r / N2);
      l = min(l, delta_r);
    } else {
      l = delta_r;
    }

    const real nut_r = c_k * l * sqrt(TKE_r);
    const real temperature_alphat_r = (one + two * l / delta_r) * nut_r;

    nut[i] = nut_r;
    temperature_alphat[i] = temperature_alphat_r;
    TKE_alphat[i] = two * nut_r;

    const real s11 = a11_r + a11_r;
    const real s22 = a22_r + a22_r;
    const real s33 = a33_r + a33_r;
    const real s12 = a12_r + a21_r;
    const real s13 = a13_r + a31_r;
    const real s23 = a23_r + a32_r;

    const real shear = nut_r * (s11 * a11_r + s12 * a12_r + s13 * a13_r
                                + s12 * a21_r + s22 * a22_r + s23 * a23_r
                                + s13 * a31_r + s23 * a32_r + s33 * a33_r);
    const real buoyancy = - (dTdx_r * g1 +
                             dTdy_r * g2 +
                             dTdz_r * g3) * temperature_alphat_r / T0;
    const real dissipation = - ((real) 0.19 + (real) 0.74 * l / delta_r)
      * sqrt(TKE_r * TKE_r * TKE_r) / l;

    TKE_source[i] = shear + buoyancy + dissipation;
  }
}

#endif // __LES_DEARDORFF_NUT_KERNEL_CL__
