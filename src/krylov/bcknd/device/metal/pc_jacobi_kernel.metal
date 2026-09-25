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
 * Metal compute kernels for the Jacobi preconditioner.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Jacobi kernel — diagonal of the stiffness matrix
 */

template <int LX>
void jacobi_impl(device float *du,
                 device const float *dxt,
                 device const float *dyt,
                 device const float *dzt,
                 device const float *G11,
                 device const float *G22,
                 device const float *G33,
                 device const float *G12,
                 device const float *G13,
                 device const float *G23,
                 const int nel,
                 uint tid, uint stride) {

  const int n = nel * LX * LX * LX;

  for (int idx = int(tid); idx < n; idx += int(stride)) {

    const int e = idx / (LX*LX*LX);
    const int ijk = idx - e*LX*LX*LX;
    const int jk = ijk / LX;
    const int i = ijk - jk * LX;
    const int k = jk / LX;
    const int j = jk - k * LX;

    float d = 0.0f;

    for (int l = 0; l < LX; l++) {
      const float g = G11[l + LX*j + LX*LX*k + LX*LX*LX*e];
      const float t = dxt[i + LX*l];
      d += g*t*t;
    }

    for (int l = 0; l < LX; l++) {
      const float g = G22[i + LX*l + LX*LX*k + LX*LX*LX*e];
      const float t = dyt[j + LX*l];
      d += g*t*t;
    }

    for (int l = 0; l < LX; l++) {
      const float g = G33[i + LX*j + LX*LX*l + LX*LX*LX*e];
      const float t = dzt[k + LX*l];
      d += g*t*t;
    }

    /* Corrections for deformed elements */
    if (i == 0 || i == LX-1) {
      d += G12[i + LX*j + LX*LX*k + LX*LX*LX*e] * dxt[i + LX*i] * dyt[j + LX*j];
      d += G13[i + LX*j + LX*LX*k + LX*LX*LX*e] * dxt[i + LX*i] * dzt[k + LX*k];
    }

    if (j == 0 || j == LX-1) {
      d += G12[i + LX*j + LX*LX*k + LX*LX*LX*e] * dyt[j + LX*j] * dxt[i + LX*i];
      d += G23[i + LX*j + LX*LX*k + LX*LX*LX*e] * dyt[j + LX*j] * dzt[k + LX*k];
    }

    if (k == 0 || k == LX-1) {
      d += G13[i + LX*j + LX*LX*k + LX*LX*LX*e] * dzt[k + LX*k] * dxt[i + LX*i];
      d += G23[i + LX*j + LX*LX*k + LX*LX*LX*e] * dzt[k + LX*k] * dyt[j + LX*j];
    }

    du[idx] = d;
  }
}

#define INSTANTIATE_JACOBI(LX)                                                  \
kernel void jacobi_kernel_lx##LX(                                               \
    device float *du[[ buffer(0) ]],                                            \
    device const float *dxt[[ buffer(1) ]],                                     \
    device const float *dyt[[ buffer(2) ]],                                     \
    device const float *dzt[[ buffer(3) ]],                                     \
    device const float *G11[[ buffer(4) ]],                                     \
    device const float *G22[[ buffer(5) ]],                                     \
    device const float *G33[[ buffer(6) ]],                                     \
    device const float *G12[[ buffer(7) ]],                                     \
    device const float *G13[[ buffer(8) ]],                                     \
    device const float *G23[[ buffer(9) ]],                                     \
    constant int &nel[[ buffer(10) ]],                                          \
    uint tid [[ thread_position_in_grid ]],                                     \
    uint stride [[ threads_per_grid ]]) {                                       \
  jacobi_impl<LX>(du, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23,             \
                  nel, tid, stride);                                           \
}

INSTANTIATE_JACOBI(2)
INSTANTIATE_JACOBI(3)
INSTANTIATE_JACOBI(4)
INSTANTIATE_JACOBI(5)
INSTANTIATE_JACOBI(6)
INSTANTIATE_JACOBI(7)
INSTANTIATE_JACOBI(8)
INSTANTIATE_JACOBI(9)
INSTANTIATE_JACOBI(10)
INSTANTIATE_JACOBI(11)
INSTANTIATE_JACOBI(12)
INSTANTIATE_JACOBI(13)
INSTANTIATE_JACOBI(14)
INSTANTIATE_JACOBI(15)
INSTANTIATE_JACOBI(16)
