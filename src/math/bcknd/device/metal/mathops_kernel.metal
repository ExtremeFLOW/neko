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
 * Metal compute kernels for field operations (mathops).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Device kernel for opchsign
 */
kernel void opchsign_kernel(device float *a1[[ buffer(0) ]],
                            device float *a2[[ buffer(1) ]],
                            device float *a3[[ buffer(2) ]],
                            constant int &gdim[[ buffer(3) ]],
                            constant int &n[[ buffer(4) ]],
                            uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  a1[idx] = -a1[idx];
  a2[idx] = -a2[idx];
  if (gdim == 3) {
    a3[idx] = -a3[idx];
  }
}

/**
 * Device kernel for opcolv
 */
kernel void opcolv_kernel(device float *a1[[ buffer(0) ]],
                          device float *a2[[ buffer(1) ]],
                          device float *a3[[ buffer(2) ]],
                          device const float *c[[ buffer(3) ]],
                          constant int &gdim[[ buffer(4) ]],
                          constant int &n[[ buffer(5) ]],
                          uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  a1[idx] = a1[idx] * c[idx];
  a2[idx] = a2[idx] * c[idx];
  if (gdim == 3) {
    a3[idx] = a3[idx] * c[idx];
  }
}

/**
 * Device kernel for opcolv3c
 */
kernel void opcolv3c_kernel(device float *a1[[ buffer(0) ]],
                            device float *a2[[ buffer(1) ]],
                            device float *a3[[ buffer(2) ]],
                            device const float *b1[[ buffer(3) ]],
                            device const float *b2[[ buffer(4) ]],
                            device const float *b3[[ buffer(5) ]],
                            device const float *c[[ buffer(6) ]],
                            constant float &d[[ buffer(7) ]],
                            constant int &gdim[[ buffer(8) ]],
                            constant int &n[[ buffer(9) ]],
                            uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  a1[idx] = b1[idx] * c[idx] * d;
  a2[idx] = b2[idx] * c[idx] * d;
  if (gdim == 3) {
    a3[idx] = b3[idx] * c[idx] * d;
  }
}

/**
 * Device kernel for opadd2cm
 */
kernel void opadd2cm_kernel(device float *a1[[ buffer(0) ]],
                            device float *a2[[ buffer(1) ]],
                            device float *a3[[ buffer(2) ]],
                            device const float *b1[[ buffer(3) ]],
                            device const float *b2[[ buffer(4) ]],
                            device const float *b3[[ buffer(5) ]],
                            constant float &c[[ buffer(6) ]],
                            constant int &gdim[[ buffer(7) ]],
                            constant int &n[[ buffer(8) ]],
                            uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  a1[idx] = a1[idx] + b1[idx] * c;
  a2[idx] = a2[idx] + b2[idx] * c;
  if (gdim == 3) {
    a3[idx] = a3[idx] + b3[idx] * c;
  }
}

/**
 * Device kernel for opadd2col
 */
kernel void opadd2col_kernel(device float *a1[[ buffer(0) ]],
                             device float *a2[[ buffer(1) ]],
                             device float *a3[[ buffer(2) ]],
                             device const float *b1[[ buffer(3) ]],
                             device const float *b2[[ buffer(4) ]],
                             device const float *b3[[ buffer(5) ]],
                             device const float *c[[ buffer(6) ]],
                             constant int &gdim[[ buffer(7) ]],
                             constant int &n[[ buffer(8) ]],
                             uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  a1[idx] = a1[idx] + b1[idx] * c[idx];
  a2[idx] = a2[idx] + b2[idx] * c[idx];
  if (gdim == 3) {
    a3[idx] = a3[idx] + b3[idx] * c[idx];
  }
}
