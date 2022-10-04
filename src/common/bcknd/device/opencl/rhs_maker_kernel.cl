/*
 Copyright (c) 2022, The Neko Authors
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

__kernel void sumab_kernel(__global real * __restrict__ u,
                           __global real * __restrict__ v,
                           __global real * __restrict__ w,
                           __global const real * __restrict__ uu,
                           __global const real * __restrict__ vv,
                           __global const real * __restrict__ ww,
                           __global const real * __restrict__ ulag1,
                           __global const real * __restrict__ ulag2,
                           __global const real * __restrict__ vlag1,
                           __global const real * __restrict__ vlag2,
                           __global const real * __restrict__ wlag1,
                           __global const real * __restrict__ wlag2,
                           const real ab1,
                           const real ab2,
                           const real ab3,
                           const int nab,
                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    u[i] = ab1 * uu[i] + ab2 * ulag1[i];
    v[i] = ab1 * vv[i] + ab2 * vlag1[i];
    w[i] = ab1 * ww[i] + ab2 * wlag1[i];
  }

  if (nab == 3) {
    for (int i = idx; i < n; i += str) {
      u[i] = u[i] + ab3 * ulag2[i];
      v[i] = v[i] + ab3 * vlag2[i];
      w[i] = w[i] + ab3 * wlag2[i];
    } 
  }
  
}

__kernel void makeext_kernel(__global real * __restrict__ abx1,
                             __global real * __restrict__ aby1,
                             __global real * __restrict__ abz1,
                             __global real * __restrict__ abx2,
                             __global real * __restrict__ aby2,
                             __global real * __restrict__ abz2,
                             __global real * __restrict__ bfx,
                             __global real * __restrict__ bfy,
                             __global real * __restrict__ bfz,
                             const real rho,
                             const real ab1,
                             const real ab2,
                             const real ab3,
                             const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    real ta1_val = ab2 * abx1[i] + ab3 * abx2[i];
    real ta2_val = ab2 * aby1[i] + ab3 * aby2[i];
    real ta3_val = ab2 * abz1[i] + ab3 * abz2[i];

    abx2[i] = abx1[i];
    aby2[i] = aby1[i];
    abz2[i] = abz1[i];
    abx1[i] = bfx[i];
    aby1[i] = bfy[i];
    abz1[i] = bfz[i];

    bfx[i] = (ab1 * bfx[i] + ta1_val) * rho;
    bfy[i] = (ab1 * bfy[i] + ta2_val) * rho;
    bfz[i] = (ab1 * bfz[i] + ta3_val) * rho;
  }
    
}

__kernel void makebdf_kernel(__global const real * __restrict__ ulag1,
                             __global const real * __restrict__ ulag2,
                             __global const real * __restrict__ vlag1,
                             __global const real * __restrict__ vlag2,
                             __global const real * __restrict__ wlag1,
                             __global const real * __restrict__ wlag2,
                             __global real * __restrict__ bfx,
                             __global real * __restrict__ bfy,
                             __global real * __restrict__ bfz,
                             __global const real * __restrict__ u,
                             __global const real * __restrict__ v,
                             __global const real * __restrict__ w,
                             __global const real * __restrict__ B,
                             const real rho,
                             const real dt,
                             const real bd2,
                             const real bd3,
                             const real bd4,
                             const int nbd,
                             const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    real tb1_val = u[i] * B[i] * bd2;
    real tb2_val = v[i] * B[i] * bd2;
    real tb3_val = w[i] * B[i] * bd2;

    real ta1_val = ulag1[i] * B[i] * bd3;
    real ta2_val = vlag1[i] * B[i] * bd3;
    real ta3_val = wlag1[i] * B[i] * bd3;

    tb1_val += ta1_val;
    tb2_val += ta2_val;
    tb3_val += ta3_val;

    if (nbd == 3) {
      tb1_val += ulag2[i] * B[i] * bd4;
      tb2_val += vlag2[i] * B[i] * bd4;
      tb3_val += wlag2[i] * B[i] * bd4;
    }
    
    bfx[i] = bfx[i] + tb1_val * (rho / dt);
    bfy[i] = bfy[i] + tb2_val * (rho / dt);
    bfz[i] = bfz[i] + tb3_val * (rho / dt);
  }
  
}


