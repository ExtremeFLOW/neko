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
 * Metal compute kernel for the cyclic (rotational) boundary rotation of
 * a velocity field at the cyclic mask points.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 * @note x, y, z and vz are bound for ABI compatibility with the other
 *       backends but are not used by this rotation.
 */

#include <metal_stdlib>
using namespace metal;

kernel void rotate_cyc_kernel(device float *vx[[ buffer(0) ]],
                              device float *vy[[ buffer(1) ]],
                              device float *vz[[ buffer(2) ]],
                              device const float *x[[ buffer(3) ]],
                              device const float *y[[ buffer(4) ]],
                              device const float *z[[ buffer(5) ]],
                              device const int *cyc_msk[[ buffer(6) ]],
                              device const float *R[[ buffer(7) ]],
                              constant int &ncyc[[ buffer(8) ]],
                              constant int &idir[[ buffer(9) ]],
                              uint tid [[ thread_position_in_grid ]]) {
  const int i = (int)tid;
  if (i >= ncyc) return;

  const int j = cyc_msk[i + 1] - 1;
  const float vxj = vx[j];
  const float vyj = vy[j];
  const float vzj = vz[j];
  device const float *Ri = R + 9 * i;

  if (idir == 1) {
    vx[j] = Ri[0] * vxj + Ri[1] * vyj + Ri[2] * vzj;
    vy[j] = Ri[3] * vxj + Ri[4] * vyj + Ri[5] * vzj;
    vz[j] = Ri[6] * vxj + Ri[7] * vyj + Ri[8] * vzj;
  } else if (idir == 0) {
    vx[j] = Ri[0] * vxj + Ri[3] * vyj + Ri[6] * vzj;
    vy[j] = Ri[1] * vxj + Ri[4] * vyj + Ri[7] * vzj;
    vz[j] = Ri[2] * vxj + Ri[5] * vyj + Ri[8] * vzj;
  }
}
