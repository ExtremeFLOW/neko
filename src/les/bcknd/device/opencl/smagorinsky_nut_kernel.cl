#ifndef __LES_SMAGORINSKY_NUT_KERNEL_CL__
#define __LES_SMAGORINSKY_NUT_KERNEL_CL__
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
 * Device kernel for smagorinsky_nut_compute
 */
__kernel void smagorinsky_nut_compute_kernel(__global const real * __restrict__ s11,
                                             __global const real * __restrict__ s22,
                                             __global const real * __restrict__ s33,
                                             __global const real * __restrict__ s12,
                                             __global const real * __restrict__ s13,
                                             __global const real * __restrict__ s23,
                                             __global const real * __restrict__ delta,
                                             __global real * __restrict__ nut,
                                             __global const real * __restrict__ mult,
                                             const real c_s,
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
    const real delta_r = delta[i];

    const real s_abs = sqrt((real) 2.0 * (s11_r * s11_r +
                                          s22_r * s22_r +
                                          s33_r * s33_r) +
                            (real) 4.0 * (s12_r * s12_r +
                                          s13_r * s13_r +
                                          s23_r * s23_r));

    nut[i] = c_s * c_s * delta_r * delta_r * s_abs * mult[i];
  }
}

#endif // __LES_SMAGORINSKY_NUT_KERNEL_CL__
