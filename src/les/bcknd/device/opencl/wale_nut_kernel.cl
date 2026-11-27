#ifndef __LES_WALE_NUT_KERNEL_CL__
#define __LES_WALE_NUT_KERNEL_CL__
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
 * Device kernel for wale_nut_compute
 */
__kernel void wale_nut_compute_kernel(__global const real * __restrict__ g11,
                                      __global const real * __restrict__ g12,
                                      __global const real * __restrict__ g13,
                                      __global const real * __restrict__ g21,
                                      __global const real * __restrict__ g22,
                                      __global const real * __restrict__ g23,
                                      __global const real * __restrict__ g31,
                                      __global const real * __restrict__ g32,
                                      __global const real * __restrict__ g33,
                                      __global const real * __restrict__ delta,
                                      __global real * __restrict__ nut,
                                      __global const real * __restrict__ mult,
                                      const real c,
                                      const real eps,
                                      const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  const real one_half = (real) 0.5;
  const real two = (real) 2.0;
  const real third = (real) 1.0 / (real) 3.0;

  for (int i = idx; i < n; i += str) {

    const real g11_r = g11[i];
    const real g12_r = g12[i];
    const real g13_r = g13[i];
    const real g21_r = g21[i];
    const real g22_r = g22[i];
    const real g23_r = g23[i];
    const real g31_r = g31[i];
    const real g32_r = g32[i];
    const real g33_r = g33[i];
    const real delta_r = delta[i];

    const real s11 = g11_r;
    const real s22 = g22_r;
    const real s33 = g33_r;
    const real s12 = one_half * (g12_r + g21_r);
    const real s13 = one_half * (g13_r + g31_r);
    const real s23 = one_half * (g23_r + g32_r);

    const real gsqr_11 = g11_r * g11_r + g12_r * g21_r + g13_r * g31_r;
    const real gsqr_12 = g11_r * g12_r + g12_r * g22_r + g13_r * g32_r;
    const real gsqr_13 = g11_r * g13_r + g12_r * g23_r + g13_r * g33_r;
    const real gsqr_21 = g21_r * g11_r + g22_r * g21_r + g23_r * g31_r;
    const real gsqr_22 = g21_r * g12_r + g22_r * g22_r + g23_r * g32_r;
    const real gsqr_23 = g21_r * g13_r + g22_r * g23_r + g23_r * g33_r;
    const real gsqr_31 = g31_r * g11_r + g32_r * g21_r + g33_r * g31_r;
    const real gsqr_32 = g31_r * g12_r + g32_r * g22_r + g33_r * g32_r;
    const real gsqr_33 = g31_r * g13_r + g32_r * g23_r + g33_r * g33_r;

    const real trace = (gsqr_11 + gsqr_22 + gsqr_33) * third;

    const real sd11 = gsqr_11 - trace;
    const real sd22 = gsqr_22 - trace;
    const real sd33 = gsqr_33 - trace;
    const real sd12 = one_half * (gsqr_12 + gsqr_21);
    const real sd13 = one_half * (gsqr_13 + gsqr_31);
    const real sd23 = one_half * (gsqr_23 + gsqr_32);

    const real Sdij_Sdij = sd11 * sd11 + sd22 * sd22 + sd33 * sd33 +
      two * (sd12 * sd12 + sd13 * sd13 + sd23 * sd23);
    const real Sij_Sij = s11 * s11 + s22 * s22 + s33 * s33 +
      two * (s12 * s12 + s13 * s13 + s23 * s23);

    const real OP_wale = sqrt(Sdij_Sdij * Sdij_Sdij * Sdij_Sdij) /
      max((sqrt(Sij_Sij * Sij_Sij * Sij_Sij * Sij_Sij * Sij_Sij)
           + sqrt(sqrt(Sdij_Sdij * Sdij_Sdij * Sdij_Sdij
                       * Sdij_Sdij * Sdij_Sdij))), eps);

    nut[i] = c * c * delta_r * delta_r * OP_wale * mult[i];
  }
}

#endif // __LES_WALE_NUT_KERNEL_CL__
