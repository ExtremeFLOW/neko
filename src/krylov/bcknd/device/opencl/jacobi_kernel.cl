#ifndef __KRYLOV_JACOBI_KERNEL_CL__
#define __KRYLOV_JACOBI_KERNEL_CL__
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

/**
 * Device kernel for jacobi
 */

#define DEFINE_JACOBI_KERNEL(LX)                                               \
__kernel void jacobi_kernel_lx##LX(__global real * __restrict__ du,            \
                                   __global const real * __restrict__ dxt,     \
                                   __global const real * __restrict__ dyt,     \
                                   __global const real * __restrict__ dzt,     \
                                   __global const real * __restrict__ G11,     \
                                   __global const real * __restrict__ G22,     \
                                   __global const real * __restrict__ G33,     \
                                   __global const real * __restrict__ G12,     \
                                   __global const real * __restrict__ G13,     \
                                   __global const real * __restrict__ G23,     \
                                   const int nel) {                            \
                                                                               \
  const int idx = get_global_id(0);                                            \
  const int e = idx / (LX*LX*LX);                                              \
  const int ijk = idx - e*LX*LX*LX;                                            \
  const int jk = ijk / LX;                                                     \
  const int i = ijk - jk * LX;                                                 \
  const int k = jk / LX;                                                       \
  const int j = jk - k * LX;                                                   \
                                                                               \
  if (e >= nel)                                                                \
    return;                                                                    \
                                                                               \
  real d = 0.0;                                                                \
                                                                               \
  for (int l = 0; l < LX; l++) {                                               \
    real g = G11[l + LX*j + LX*LX*k + LX*LX*LX*e];                             \
    real t = dxt[i + LX*l];                                                    \
    d += g*t*t;                                                                \
  }                                                                            \
                                                                               \
  for (int l = 0; l < LX; l++) {                                               \
    real g = G22[i + LX*l + LX*LX*k + LX*LX*LX*e];                             \
    real t = dyt[j + LX*l];                                                    \
    d += g*t*t;                                                                \
  }                                                                            \
                                                                               \
  for (int l = 0; l < LX; l++) {                                               \
    real g = G33[i + LX*j + LX*LX*l + LX*LX*LX*e];                             \
    real t = dzt[k + LX*l];                                                    \
    d += g*t*t;                                                                \
  }                                                                            \
                                                                               \
  /* Corrections for deformed elements */                                      \
  if (i == 0 || i == LX-1) {                                                   \
    d += G12[i + LX*j + LX*LX*k + LX*LX*LX*e] * dxt[i + LX*i] * dyt[j + LX*j]; \
    d += G13[i + LX*j + LX*LX*k + LX*LX*LX*e] * dxt[i + LX*i] * dzt[k + LX*k]; \
  }                                                                            \
                                                                               \
  if (j == 0 || j == LX-1) {                                                   \
    d += G12[i + LX*j + LX*LX*k + LX*LX*LX*e] * dyt[j + LX*j] * dxt[i + LX*i]; \
    d += G23[i + LX*j + LX*LX*k + LX*LX*LX*e] * dyt[j + LX*j] * dzt[k + LX*k]; \
  }                                                                            \
                                                                               \
  if (k == 0 || k == LX-1) {                                                   \
    d += G13[i + LX*j + LX*LX*k + LX*LX*LX*e] * dzt[k + LX*k] * dxt[i + LX*i]; \
    d += G23[i + LX*j + LX*LX*k + LX*LX*LX*e] * dzt[k + LX*k] * dyt[j + LX*j]; \
  }                                                                            \
                                                                               \
  du[idx] = d;                                                                 \
}

DEFINE_JACOBI_KERNEL(2)
DEFINE_JACOBI_KERNEL(3)
DEFINE_JACOBI_KERNEL(4)
DEFINE_JACOBI_KERNEL(5)
DEFINE_JACOBI_KERNEL(6)
DEFINE_JACOBI_KERNEL(7)
DEFINE_JACOBI_KERNEL(8)
DEFINE_JACOBI_KERNEL(9)
DEFINE_JACOBI_KERNEL(10)
DEFINE_JACOBI_KERNEL(11)
DEFINE_JACOBI_KERNEL(12)
DEFINE_JACOBI_KERNEL(13)
DEFINE_JACOBI_KERNEL(14)
DEFINE_JACOBI_KERNEL(15)
DEFINE_JACOBI_KERNEL(16)


#endif // __KRYLOV_JACOBI_KERNEL_CL__
