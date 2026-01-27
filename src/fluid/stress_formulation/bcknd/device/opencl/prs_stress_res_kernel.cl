/*
 Copyright (c) 2025, The Neko Authors
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

__kernel
void prs_stress_res_part1_kernel(__global real * __restrict__ ta1,
                                 __global real * __restrict__ ta2,
                                 __global real * __restrict__ ta3,
                                 __global real * __restrict__ wa1,
                                 __global real * __restrict__ wa2,
                                 __global real * __restrict__ wa3,
                                 __global const real * __restrict__ s11,
                                 __global const real * __restrict__ s22,
                                 __global const real * __restrict__ s33,
                                 __global const real * __restrict__ s12,
                                 __global const real * __restrict__ s13,
                                 __global const real * __restrict__ s23,
                                 __global const real * __restrict__ f_u,
                                 __global const real * __restrict__ f_v,
                                 __global const real * __restrict__ f_w,
                                 __global const real * __restrict__ B,
                                 __global const real * __restrict__ rho,
                                 const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    wa1[i] -= 2.0 * (ta1[i] * s11[i]
                   + ta2[i] * s12[i]
                   + ta3[i] * s13[i]);
    wa2[i] -= 2.0 * (ta1[i] * s12[i]
                   + ta2[i] * s22[i]
                   + ta3[i] * s23[i]);
    wa3[i] -= 2.0 * (ta1[i] * s13[i]
                   + ta2[i] * s23[i]
                   + ta3[i] * s33[i]);

    ta1[i] = (f_u[i] / rho[i]) - ((wa1[i] / rho[i]) * B[i]);
    ta2[i] = (f_v[i] / rho[i]) - ((wa2[i] / rho[i]) * B[i]);
    ta3[i] = (f_w[i] / rho[i]) - ((wa3[i] / rho[i]) * B[i]);
  }

}

__kernel
void prs_stress_res_part3_kernel(__global real * __restrict__ p_res,
                                 __global const real * __restrict__ ta1,
                                 __global const real * __restrict__ ta2,
                                 __global const real * __restrict__ ta3,
                                 __global const real * __restrict__ wa1,
                                 __global const real * __restrict__ wa2,
                                 __global const real * __restrict__ wa3,
                                 const real dtbd,
                                 const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    p_res[i] = p_res[i] - (dtbd * (ta1[i] + ta2[i] + ta3[i]))
      - (wa1[i] + wa2[i] + wa3[i]);
  }
}
