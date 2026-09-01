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
 * Metal compute kernels for the rhs makers (sumab/ext/bdf/oifs).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Sum the Adams-Bashforth extrapolation of the velocity field.
 */
kernel void sumab_kernel(device float *u[[ buffer(0) ]],
                         device float *v[[ buffer(1) ]],
                         device float *w[[ buffer(2) ]],
                         device const float *uu[[ buffer(3) ]],
                         device const float *vv[[ buffer(4) ]],
                         device const float *ww[[ buffer(5) ]],
                         device const float *ulag1[[ buffer(6) ]],
                         device const float *ulag2[[ buffer(7) ]],
                         device const float *vlag1[[ buffer(8) ]],
                         device const float *vlag2[[ buffer(9) ]],
                         device const float *wlag1[[ buffer(10) ]],
                         device const float *wlag2[[ buffer(11) ]],
                         constant float &ab1[[ buffer(12) ]],
                         constant float &ab2[[ buffer(13) ]],
                         constant float &ab3[[ buffer(14) ]],
                         constant int &nab[[ buffer(15) ]],
                         constant int &n[[ buffer(16) ]],
                         uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  u[idx] = ab1 * uu[idx] + ab2 * ulag1[idx];
  v[idx] = ab1 * vv[idx] + ab2 * vlag1[idx];
  w[idx] = ab1 * ww[idx] + ab2 * wlag1[idx];
  if (nab == 3) {
    u[idx] = u[idx] + ab3 * ulag2[idx];
    v[idx] = v[idx] + ab3 * vlag2[idx];
    w[idx] = w[idx] + ab3 * wlag2[idx];
  }
}

/**
 * Extrapolate the velocity source terms (makeext).
 */
kernel void makeext_kernel(device float *abx1[[ buffer(0) ]],
                           device float *aby1[[ buffer(1) ]],
                           device float *abz1[[ buffer(2) ]],
                           device float *abx2[[ buffer(3) ]],
                           device float *aby2[[ buffer(4) ]],
                           device float *abz2[[ buffer(5) ]],
                           device float *bfx[[ buffer(6) ]],
                           device float *bfy[[ buffer(7) ]],
                           device float *bfz[[ buffer(8) ]],
                           constant float &rho[[ buffer(9) ]],
                           constant float &ab1[[ buffer(10) ]],
                           constant float &ab2[[ buffer(11) ]],
                           constant float &ab3[[ buffer(12) ]],
                           constant int &n[[ buffer(13) ]],
                           uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float ta1_val = ab2 * abx1[idx] + ab3 * abx2[idx];
  const float ta2_val = ab2 * aby1[idx] + ab3 * aby2[idx];
  const float ta3_val = ab2 * abz1[idx] + ab3 * abz2[idx];

  abx2[idx] = abx1[idx];
  aby2[idx] = aby1[idx];
  abz2[idx] = abz1[idx];
  abx1[idx] = bfx[idx];
  aby1[idx] = bfy[idx];
  abz1[idx] = bfz[idx];

  bfx[idx] = (ab1 * bfx[idx] + ta1_val) * rho;
  bfy[idx] = (ab1 * bfy[idx] + ta2_val) * rho;
  bfz[idx] = (ab1 * bfz[idx] + ta3_val) * rho;
}

/**
 * Extrapolate the scalar source term (scalar_makeext).
 */
kernel void scalar_makeext_kernel(device float *fs_lag[[ buffer(0) ]],
                                  device float *fs_laglag[[ buffer(1) ]],
                                  device float *fs[[ buffer(2) ]],
                                  constant float &ext1[[ buffer(3) ]],
                                  constant float &ext2[[ buffer(4) ]],
                                  constant float &ext3[[ buffer(5) ]],
                                  constant int &n[[ buffer(6) ]],
                                  uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float ta1_val = ext2 * fs_lag[idx] + ext3 * fs_laglag[idx];
  fs_laglag[idx] = fs_lag[idx];
  fs_lag[idx] = fs[idx];
  fs[idx] = ext1 * fs[idx] + ta1_val;
}

/**
 * BDF contribution to the velocity source terms (makebdf).
 */
kernel void makebdf_kernel(device const float *ulag1[[ buffer(0) ]],
                           device const float *ulag2[[ buffer(1) ]],
                           device const float *vlag1[[ buffer(2) ]],
                           device const float *vlag2[[ buffer(3) ]],
                           device const float *wlag1[[ buffer(4) ]],
                           device const float *wlag2[[ buffer(5) ]],
                           device float *bfx[[ buffer(6) ]],
                           device float *bfy[[ buffer(7) ]],
                           device float *bfz[[ buffer(8) ]],
                           device const float *u[[ buffer(9) ]],
                           device const float *v[[ buffer(10) ]],
                           device const float *w[[ buffer(11) ]],
                           device const float *B[[ buffer(12) ]],
                           constant float &rho[[ buffer(13) ]],
                           constant float &dt[[ buffer(14) ]],
                           constant float &bd2[[ buffer(15) ]],
                           constant float &bd3[[ buffer(16) ]],
                           constant float &bd4[[ buffer(17) ]],
                           constant int &nbd[[ buffer(18) ]],
                           constant int &n[[ buffer(19) ]],
                           uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  float tb1_val = u[idx] * B[idx] * bd2;
  float tb2_val = v[idx] * B[idx] * bd2;
  float tb3_val = w[idx] * B[idx] * bd2;

  tb1_val += ulag1[idx] * B[idx] * bd3;
  tb2_val += vlag1[idx] * B[idx] * bd3;
  tb3_val += wlag1[idx] * B[idx] * bd3;

  if (nbd == 3) {
    tb1_val += ulag2[idx] * B[idx] * bd4;
    tb2_val += vlag2[idx] * B[idx] * bd4;
    tb3_val += wlag2[idx] * B[idx] * bd4;
  }

  bfx[idx] = bfx[idx] + tb1_val * (rho / dt);
  bfy[idx] = bfy[idx] + tb2_val * (rho / dt);
  bfz[idx] = bfz[idx] + tb3_val * (rho / dt);
}

/**
 * BDF contribution to the scalar source term (scalar_makebdf).
 */
kernel void scalar_makebdf_kernel(device const float *s_lag[[ buffer(0) ]],
                                  device const float *s_laglag[[ buffer(1) ]],
                                  device float *fs[[ buffer(2) ]],
                                  device const float *s[[ buffer(3) ]],
                                  device const float *B[[ buffer(4) ]],
                                  device const float *rho_cp[[ buffer(5) ]],
                                  constant float &dt[[ buffer(6) ]],
                                  constant float &bd2[[ buffer(7) ]],
                                  constant float &bd3[[ buffer(8) ]],
                                  constant float &bd4[[ buffer(9) ]],
                                  constant int &nbd[[ buffer(10) ]],
                                  constant int &n[[ buffer(11) ]],
                                  uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  float tb1_val = s[idx] * B[idx] * bd2;
  tb1_val += s_lag[idx] * B[idx] * bd3;
  if (nbd == 3) {
    tb1_val += s_laglag[idx] * B[idx] * bd4;
  }
  fs[idx] = fs[idx] + tb1_val * (rho_cp[idx] / dt);
}

/**
 * OIFS contribution to the velocity source terms (makeoifs).
 */
kernel void makeoifs_kernel(device const float *phi_x[[ buffer(0) ]],
                            device const float *phi_y[[ buffer(1) ]],
                            device const float *phi_z[[ buffer(2) ]],
                            device float *bf_x[[ buffer(3) ]],
                            device float *bf_y[[ buffer(4) ]],
                            device float *bf_z[[ buffer(5) ]],
                            constant float &rho[[ buffer(6) ]],
                            constant float &dt[[ buffer(7) ]],
                            constant int &n[[ buffer(8) ]],
                            uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  bf_x[idx] = bf_x[idx] + phi_x[idx] * (rho / dt);
  bf_y[idx] = bf_y[idx] + phi_y[idx] * (rho / dt);
  bf_z[idx] = bf_z[idx] + phi_z[idx] * (rho / dt);
}

/**
 * OIFS contribution to the scalar source term (scalar_makeoifs).
 */
kernel void scalar_makeoifs_kernel(device const float *phi_s[[ buffer(0) ]],
                                   device float *bf_s[[ buffer(1) ]],
                                   device const float *rho_cp[[ buffer(2) ]],
                                   constant float &dt[[ buffer(3) ]],
                                   constant int &n[[ buffer(4) ]],
                                   uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  bf_s[idx] = bf_s[idx] + phi_s[idx] * (rho_cp[idx] / dt);
}
