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
 * Metal compute kernels for the lambda2 vortex criterion.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Median eigenvalue of the symmetric tensor S^2 + Omega^2
 */

/**
 * Return the middle (lambda2) eigenvalue of S^2 + Omega^2, where S and
 * Omega are the symmetric and antisymmetric parts of the velocity
 * gradient tensor.
 *
 * @note In FP32 the closed-form cubic-root expression is fragile: a
 * zero/near-zero strain field gives 0/0, and rounding can push the
 * acos() argument just outside [-1, 1]. Under Metal's fast-math these
 * do not reliably yield NaN (the CUDA path relies on NaN comparisons
 * falling through to 0), so we guard the degenerate case explicitly
 * and clamp the acos() argument.
 */
static inline float eigen_val_calc(float grad11, float grad12, float grad13,
                                   float grad21, float grad22, float grad23,
                                   float grad31, float grad32, float grad33) {
  float s11 = grad11;
  float s22 = grad22;
  float s33 = grad33;
  float s12 = 0.5f * (grad12 + grad21);
  float s13 = 0.5f * (grad13 + grad31);
  float s23 = 0.5f * (grad23 + grad32);

  float o12 = 0.5f * (grad12 - grad21);
  float o13 = 0.5f * (grad13 - grad31);
  float o23 = 0.5f * (grad23 - grad32);

  float a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13;
  float a12 = s11*s12 + s12*s22 + s13*s23 - o13*o23;
  float a13 = s11*s13 + s12*s23 + s13*s33 + o12*o23;

  float a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23;
  float a23 = s12*s13 + s22*s23 + s23*s33 - o12*o13;
  float a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23;

  float B = -(a11 + a22 + a33);
  float C = -(a12*a12 + a13*a13 + a23*a23 - a11*a22 - a11*a33 - a22*a33);
  float D = -(2.0f * a12 * a13 * a23 - a11 * a23*a23 - a22 * a13*a13
              - a33 * a12*a12 + a11 * a22 * a33);

  float q = (3.0f * C - B*B) / 9.0f;
  float r = (9.0f * C * B - 27.0f * D - 2.0f * B*B*B) / 54.0f;

  // Degenerate strain (q >= 0, including the zero-gradient 0/0 case):
  // no well-defined vortex, return 0 as the CUDA/CPU paths effectively do.
  float denom = sqrt(-q*q*q);
  if (!(denom > 0.0f)) {
    return 0.0f;
  }

  float theta = acos(clamp(r / denom, -1.0f, 1.0f));
  float pi = M_PI_F;
  float two_sqrt_negq = 2.0f * sqrt(-q);
  float eigen1 = two_sqrt_negq * cos(theta / 3.0f) - B / 3.0f;
  float eigen2 = two_sqrt_negq * cos((theta + 2.0f * pi) / 3.0f) - B / 3.0f;
  float eigen3 = two_sqrt_negq * cos((theta + 4.0f * pi) / 3.0f) - B / 3.0f;

  if (eigen1 <= eigen2 && eigen2 <= eigen3)
    return eigen2;
  else if (eigen3 <= eigen2 && eigen2 <= eigen1)
    return eigen2;
  else if (eigen1 <= eigen3 && eigen3 <= eigen2)
    return eigen3;
  else if (eigen2 <= eigen3 && eigen3 <= eigen1)
    return eigen3;
  else if (eigen2 <= eigen1 && eigen1 <= eigen3)
    return eigen1;
  else if (eigen3 <= eigen1 && eigen1 <= eigen2)
    return eigen1;
  return 0.0f;
}

/*
 * lambda2 — kstep formulation, one element per (LX x LX) threadgroup
 */

template <int LX>
void lambda2_kstep_impl(device float *lambda2,
                        device const float *u,
                        device const float *v,
                        device const float *w,
                        device const float *dx,
                        device const float *dy,
                        device const float *dz,
                        device const float *drdx,
                        device const float *dsdx,
                        device const float *dtdx,
                        device const float *drdy,
                        device const float *dsdy,
                        device const float *dtdy,
                        device const float *drdz,
                        device const float *dsdz,
                        device const float *dtdz,
                        device const float *jacinv,
                        threadgroup float *shu,
                        threadgroup float *shv,
                        threadgroup float *shw,
                        threadgroup float *shdx,
                        threadgroup float *shdy,
                        threadgroup float *shdz,
                        uint e, uint i, uint j) {

  const int ij = i + j * LX;
  const int ele = e * LX * LX * LX;

  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];

  float ru[LX];
  float rv[LX];
  float rw[LX];
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k * LX * LX + ele];
    rv[k] = v[ij + k * LX * LX + ele];
    rw[k] = w[ij + k * LX * LX + ele];
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k * LX * LX;
    const float jinv = jacinv[ijk + ele];

    float ttmpu = 0.0f;
    float ttmpv = 0.0f;
    float ttmpw = 0.0f;
    shu[ij] = ru[k];
    shv[ij] = rv[k];
    shw[ij] = rw[k];
    for (int l = 0; l < LX; l++) {
      ttmpu += shdz[k + l * LX] * ru[l];
      ttmpv += shdz[k + l * LX] * rv[l];
      ttmpw += shdz[k + l * LX] * rw[l];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float rtmpu = 0.0f;
    float stmpu = 0.0f;
    float rtmpv = 0.0f;
    float stmpv = 0.0f;
    float rtmpw = 0.0f;
    float stmpw = 0.0f;
    for (int l = 0; l < LX; l++) {
      rtmpu += shdx[i + l * LX] * shu[l + j * LX];
      stmpu += shdy[j + l * LX] * shu[i + l * LX];
      rtmpv += shdx[i + l * LX] * shv[l + j * LX];
      stmpv += shdy[j + l * LX] * shv[i + l * LX];
      rtmpw += shdx[i + l * LX] * shw[l + j * LX];
      stmpw += shdy[j + l * LX] * shw[i + l * LX];
    }

    float grad11 = jinv * (drdx[ijk + ele] * rtmpu
                           + dsdx[ijk + ele] * stmpu
                           + dtdx[ijk + ele] * ttmpu);
    float grad12 = jinv * (drdy[ijk + ele] * rtmpu
                           + dsdy[ijk + ele] * stmpu
                           + dtdy[ijk + ele] * ttmpu);
    float grad13 = jinv * (drdz[ijk + ele] * rtmpu
                           + dsdz[ijk + ele] * stmpu
                           + dtdz[ijk + ele] * ttmpu);
    float grad21 = jinv * (drdx[ijk + ele] * rtmpv
                           + dsdx[ijk + ele] * stmpv
                           + dtdx[ijk + ele] * ttmpv);
    float grad22 = jinv * (drdy[ijk + ele] * rtmpv
                           + dsdy[ijk + ele] * stmpv
                           + dtdy[ijk + ele] * ttmpv);
    float grad23 = jinv * (drdz[ijk + ele] * rtmpv
                           + dsdz[ijk + ele] * stmpv
                           + dtdz[ijk + ele] * ttmpv);
    float grad31 = jinv * (drdx[ijk + ele] * rtmpw
                           + dsdx[ijk + ele] * stmpw
                           + dtdx[ijk + ele] * ttmpw);
    float grad32 = jinv * (drdy[ijk + ele] * rtmpw
                           + dsdy[ijk + ele] * stmpw
                           + dtdy[ijk + ele] * ttmpw);
    float grad33 = jinv * (drdz[ijk + ele] * rtmpw
                           + dsdz[ijk + ele] * stmpw
                           + dtdz[ijk + ele] * ttmpw);

    lambda2[ijk + ele] = eigen_val_calc(grad11, grad12, grad13,
                                        grad21, grad22, grad23,
                                        grad31, grad32, grad33);
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }
}

#define INSTANTIATE_LAMBDA2(LX)                                                \
kernel void lambda2_kernel_lx##LX(                                            \
    device float *lambda2[[ buffer(0) ]],                                     \
    device const float *u[[ buffer(1) ]],                                     \
    device const float *v[[ buffer(2) ]],                                     \
    device const float *w[[ buffer(3) ]],                                     \
    device const float *dx[[ buffer(4) ]],                                    \
    device const float *dy[[ buffer(5) ]],                                    \
    device const float *dz[[ buffer(6) ]],                                    \
    device const float *drdx[[ buffer(7) ]],                                  \
    device const float *dsdx[[ buffer(8) ]],                                  \
    device const float *dtdx[[ buffer(9) ]],                                  \
    device const float *drdy[[ buffer(10) ]],                                 \
    device const float *dsdy[[ buffer(11) ]],                                 \
    device const float *dtdy[[ buffer(12) ]],                                 \
    device const float *drdz[[ buffer(13) ]],                                 \
    device const float *dsdz[[ buffer(14) ]],                                 \
    device const float *dtdz[[ buffer(15) ]],                                 \
    device const float *jacinv[[ buffer(16) ]],                               \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                            \
    uint2 tid [[ thread_position_in_threadgroup ]]) {                         \
  threadgroup float shu[LX * LX];                                             \
  threadgroup float shv[LX * LX];                                             \
  threadgroup float shw[LX * LX];                                             \
  threadgroup float shdx[LX * LX];                                            \
  threadgroup float shdy[LX * LX];                                            \
  threadgroup float shdz[LX * LX];                                            \
  lambda2_kstep_impl<LX>(lambda2, u, v, w, dx, dy, dz,                        \
                         drdx, dsdx, dtdx, drdy, dsdy, dtdy,                  \
                         drdz, dsdz, dtdz, jacinv,                            \
                         shu, shv, shw, shdx, shdy, shdz,                     \
                         gid2.x, tid.x, tid.y);                              \
}

INSTANTIATE_LAMBDA2(2)
INSTANTIATE_LAMBDA2(3)
INSTANTIATE_LAMBDA2(4)
INSTANTIATE_LAMBDA2(5)
INSTANTIATE_LAMBDA2(6)
INSTANTIATE_LAMBDA2(7)
INSTANTIATE_LAMBDA2(8)
INSTANTIATE_LAMBDA2(9)
INSTANTIATE_LAMBDA2(10)
INSTANTIATE_LAMBDA2(11)
INSTANTIATE_LAMBDA2(12)
INSTANTIATE_LAMBDA2(13)
INSTANTIATE_LAMBDA2(14)
INSTANTIATE_LAMBDA2(15)
INSTANTIATE_LAMBDA2(16)
