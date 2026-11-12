#ifndef __MATH_ROTATE_KERNEL_CL__
#define __MATH_ROTATE_KERNEL_CL__
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
 * Device kernel for cyclic boundary condition rotation
 */
__kernel void rotate_cyc_kernel(__global real *vx,
                                __global real *vy,
                                __global real *vz,
                                __global const real *x,
                                __global const real *y,
                                __global const real *z,
                                __global const int *cyc_msk,
                                __global const real *R,
                                const int ncyc,
                                const int idir) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < ncyc; i += str) {
    const int j = cyc_msk[i + 1] - 1;
    const real vxj = vx[j];
    const real vyj = vy[j];
    const real vzj = vz[j];
    __global const real *Ri = R + 9 * i;
    if (idir == 1) {
      vx[j] = Ri[0] * vxj + Ri[1] * vyj + Ri[2] * vzj;
      vy[j] = Ri[3] * vxj + Ri[4] * vyj + Ri[5] * vzj;
      vz[j] = Ri[6] * vxj + Ri[7] * vyj + Ri[8] * vzj;
    }
    else if (idir == 0) {
      vx[j] = Ri[0] * vxj + Ri[3] * vyj + Ri[6] * vzj;
      vy[j] = Ri[1] * vxj + Ri[4] * vyj + Ri[7] * vzj;
      vz[j] = Ri[2] * vxj + Ri[5] * vyj + Ri[8] * vzj;
    }
  }
}

#endif // __MATH_ROTATE_KERNEL_CL__
