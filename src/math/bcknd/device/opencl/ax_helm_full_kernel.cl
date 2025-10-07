#ifndef __MATH_AX_HELM_FULL_KERNEL_CL__
#define __MATH_AX_HELM_FULL_KERNEL_CL__
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

/**
 * Device kernel for Ax helm full
 */

#define DEFINE_AX_HELM_FULL_KERNEL(LX)                                         \
__kernel void                                                                  \
ax_helm_stress_kernel_full_lx##LX(__global real * __restrict__ au,             \
                                  __global real * __restrict__ av,             \
                                  __global real * __restrict__ aw,             \
                                  __global const real * __restrict__ u,        \
                                  __global const real * __restrict__ v,        \
                                  __global const real * __restrict__ w,        \
                                  __global const real * __restrict__ dx,       \
                                  __global const real * __restrict__ dy,       \
                                  __global const real * __restrict__ dz,       \
                                  __global const real * __restrict__ h1,       \
                                  __global const real * __restrict__ drdx,     \
                                  __global const real * __restrict__ drdy,     \
                                  __global const real * __restrict__ drdz,     \
                                  __global const real * __restrict__ dsdx,     \
                                  __global const real * __restrict__ dsdy,     \
                                  __global const real * __restrict__ dsdz,     \
                                  __global const real * __restrict__ dtdx,     \
                                  __global const real * __restrict__ dtdy,     \
                                  __global const real * __restrict__ dtdz,     \
                                  __global const real * __restrict__ jacinv,   \
                                  __global const real * __restrict__ weight3) {\
                                                                               \
  __local real shdx[LX * LX];                                                  \
  __local real shdy[LX * LX];                                                  \
  __local real shdz[LX * LX];                                                  \
                                                                               \
  __local real shu[LX * LX];                                                   \
  __local real shur[LX * LX];                                                  \
  __local real shus[LX * LX];                                                  \
                                                                               \
  __local real shv[LX * LX];                                                   \
  __local real shvr[LX * LX];                                                  \
  __local real shvs[LX * LX];                                                  \
                                                                               \
  __local real shw[LX * LX];                                                   \
  __local real shwr[LX * LX];                                                  \
  __local real shws[LX * LX];                                                  \
                                                                               \
  real ru[LX];                                                                 \
  real rv[LX];                                                                 \
  real rw[LX];                                                                 \
                                                                               \
  real ruw[LX];                                                                \
  real rvw[LX];                                                                \
  real rww[LX];                                                                \
                                                                               \
  real rut;                                                                    \
  real rvt;                                                                    \
  real rwt;                                                                    \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int j = get_local_id(1);                                               \
  const int i = get_local_id(0);                                               \
  const int ij = i + j*LX;                                                     \
  const int ele = e*LX*LX*LX;                                                  \
                                                                               \
  shdx[ij] = dx[ij];                                                           \
  shdy[ij] = dy[ij];                                                           \
  shdz[ij] = dz[ij];                                                           \
                                                                               \
  for(int k = 0; k < LX; ++k){                                                 \
    ru[k] = u[ij + k*LX*LX + ele];                                             \
    ruw[k] = 0.0;                                                              \
                                                                               \
    rv[k] = v[ij + k*LX*LX + ele];                                             \
    rvw[k] = 0.0;                                                              \
                                                                               \
    rw[k] = w[ij + k*LX*LX + ele];                                             \
    rww[k] = 0.0;                                                              \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int k = 0; k < LX; ++k){                                                \
    const int ijk = ij + k*LX*LX;                                              \
    const real drdx_local = drdx[ijk+ele];                                     \
    const real drdy_local = drdy[ijk+ele];                                     \
    const real drdz_local = drdz[ijk+ele];                                     \
    const real dsdx_local = dsdx[ijk+ele];                                     \
    const real dsdy_local = dsdy[ijk+ele];                                     \
    const real dsdz_local = dsdz[ijk+ele];                                     \
    const real dtdx_local = dtdx[ijk+ele];                                     \
    const real dtdy_local = dtdy[ijk+ele];                                     \
    const real dtdz_local = dtdz[ijk+ele];                                     \
    const real dj  = h1[ijk+ele] *                                             \
                     weight3[ijk] *                                            \
                     jacinv[ijk+ele];                                          \
                                                                               \
    real uttmp = 0.0;                                                          \
    real vttmp = 0.0;                                                          \
    real wttmp = 0.0;                                                          \
    shu[ij] = ru[k];                                                           \
    shv[ij] = rv[k];                                                           \
    shw[ij] = rw[k];                                                           \
    for (int l = 0; l < LX; l++){                                              \
      uttmp += shdz[k+l*LX] * ru[l];                                           \
      vttmp += shdz[k+l*LX] * rv[l];                                           \
      wttmp += shdz[k+l*LX] * rw[l];                                           \
    }                                                                          \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
                                                                               \
    real urtmp = 0.0;                                                          \
    real ustmp = 0.0;                                                          \
                                                                               \
    real vrtmp = 0.0;                                                          \
    real vstmp = 0.0;                                                          \
                                                                               \
    real wrtmp = 0.0;                                                          \
    real wstmp = 0.0;                                                          \
                                                                               \
    for (int l = 0; l < LX; l++){                                              \
      urtmp += shdx[i+l*LX] * shu[l+j*LX];                                     \
      ustmp += shdy[j+l*LX] * shu[i+l*LX];                                     \
                                                                               \
      vrtmp += shdx[i+l*LX] * shv[l+j*LX];                                     \
      vstmp += shdy[j+l*LX] * shv[i+l*LX];                                     \
                                                                               \
      wrtmp += shdx[i+l*LX] * shw[l+j*LX];                                     \
      wstmp += shdy[j+l*LX] * shw[i+l*LX];                                     \
    }                                                                          \
                                                                               \
    real u1 = 0.0;                                                             \
    real u2 = 0.0;                                                             \
    real u3 = 0.0;                                                             \
    real v1 = 0.0;                                                             \
    real v2 = 0.0;                                                             \
    real v3 = 0.0;                                                             \
    real w1 = 0.0;                                                             \
    real w2 = 0.0;                                                             \
    real w3 = 0.0;                                                             \
                                                                               \
    u1 = urtmp * drdx_local +                                                  \
         ustmp * dsdx_local +                                                  \
         uttmp * dtdx_local;                                                   \
    u2 = urtmp * drdy_local +                                                  \
         ustmp * dsdy_local +                                                  \
         uttmp * dtdy_local;                                                   \
    u3 = urtmp * drdz_local +                                                  \ 
         ustmp * dsdz_local +                                                  \ 
         uttmp * dtdz_local;                                                   \
                                                                               \
    v1 = vrtmp * drdx_local +                                                  \ 
         vstmp * dsdx_local +                                                  \
         vttmp * dtdx_local;                                                   \
    v2 = vrtmp * drdy_local +                                                  \
         vstmp * dsdy_local +                                                  \ 
         vttmp * dtdy_local;                                                   \
    v3 = vrtmp * drdz_local +                                                  \
         vstmp * dsdz_local +                                                  \
         vttmp * dtdz_local;                                                   \
                                                                               \
    w1 = wrtmp * drdx_local +                                                  \
         wstmp * dsdx_local +                                                  \
         wttmp * dtdx_local;                                                   \
    w2 = wrtmp * drdy_local +                                                  \
         wstmp * dsdy_local +                                                  \
         wttmp * dtdy_local;                                                   \
    w3 = wrtmp * drdz_local +                                                  \
         wstmp * dsdz_local +                                                  \
         wttmp * dtdz_local;                                                   \
                                                                               \
    real s11 = 0.0;                                                            \
    real s12 = 0.0;                                                            \
    real s13 = 0.0;                                                            \
    real s21 = 0.0;                                                            \
    real s22 = 0.0;                                                            \
    real s23 = 0.0;                                                            \
    real s31 = 0.0;                                                            \
    real s32 = 0.0;                                                            \
    real s33 = 0.0;                                                            \
                                                                               \
    s11 = dj*(u1 + u1);                                                        \
    s12 = dj*(u2 + v1);                                                        \
    s13 = dj*(u3 + w1);                                                        \
    s21 = dj*(v1 + u2);                                                        \
    s22 = dj*(v2 + v2);                                                        \
    s23 = dj*(v3 + w2);                                                        \
    s31 = dj*(w1 + u3);                                                        \
    s32 = dj*(w2 + v3);                                                        \
    s33 = dj*(w3 + w3);                                                        \
                                                                               \
    shur[ij] = drdx_local * s11 +                                              \
               drdy_local * s12 +                                              \
               drdz_local * s13;                                               \
    shus[ij] = dsdx_local * s11 +                                              \
               dsdy_local * s12 +                                              \
               dsdz_local * s13;                                               \
    rut =      dtdx_local * s11 +                                              \
               dtdy_local * s12 +                                              \
               dtdz_local * s13;                                               \
                                                                               \
    shvr[ij] = drdx_local * s21 +                                              \
               drdy_local * s22 +                                              \
               drdz_local * s23;                                               \
    shvs[ij] = dsdx_local * s21 +                                              \
               dsdy_local * s22 +                                              \
               dsdz_local * s23;                                               \
    rvt =      dtdx_local * s21 +                                              \
               dtdy_local * s22 +                                              \
               dtdz_local * s23;                                               \
                                                                               \
    shwr[ij] = drdx_local * s31 +                                              \
               drdy_local * s32 +                                              \
               drdz_local * s33;                                               \
    shws[ij] = dsdx_local * s31 +                                              \
               dsdy_local * s32 +                                              \
               dsdz_local * s33;                                               \
    rwt =      dtdx_local * s31 +                                              \
               dtdy_local * s32 +                                              \
               dtdz_local * s33;                                               \
                                                                               \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
                                                                               \
    real uwijke = 0.0;                                                         \
    real vwijke = 0.0;                                                         \
    real wwijke = 0.0;                                                         \
                                                                               \
    for (int l = 0; l < LX; l++){                                              \
      uwijke += shur[l+j*LX] * shdx[l+i*LX];                                   \
      ruw[l] += rut * shdz[k+l*LX];                                            \
      uwijke += shus[i+l*LX] * shdy[l + j*LX];                                 \
                                                                               \
      vwijke += shvr[l+j*LX] * shdx[l+i*LX];                                   \
      rvw[l] += rvt * shdz[k+l*LX];                                            \
      vwijke += shvs[i+l*LX] * shdy[l + j*LX];                                 \
                                                                               \
      wwijke += shwr[l+j*LX] * shdx[l+i*LX];                                   \
      rww[l] += rwt * shdz[k+l*LX];                                            \
      wwijke += shws[i+l*LX] * shdy[l + j*LX];                                 \
    }                                                                          \
    ruw[k] += uwijke;                                                          \
    rvw[k] += vwijke;                                                          \
    rww[k] += wwijke;                                                          \
  }                                                                            \
                                                                               \
  for (int k = 0; k < LX; ++k){                                                \
   au[ij + k*LX*LX + ele] = ruw[k];                                            \
   av[ij + k*LX*LX + ele] = rvw[k];                                            \
   aw[ij + k*LX*LX + ele] = rww[k];                                            \
  }                                                                            \
}
DEFINE_AX_HELM_FULL_KERNEL(2)
DEFINE_AX_HELM_FULL_KERNEL(3)
DEFINE_AX_HELM_FULL_KERNEL(4)
DEFINE_AX_HELM_FULL_KERNEL(5)
DEFINE_AX_HELM_FULL_KERNEL(6)
DEFINE_AX_HELM_FULL_KERNEL(7)
DEFINE_AX_HELM_FULL_KERNEL(8)
DEFINE_AX_HELM_FULL_KERNEL(9)
DEFINE_AX_HELM_FULL_KERNEL(10)
DEFINE_AX_HELM_FULL_KERNEL(11)
DEFINE_AX_HELM_FULL_KERNEL(12)
DEFINE_AX_HELM_FULL_KERNEL(13)
DEFINE_AX_HELM_FULL_KERNEL(14)
DEFINE_AX_HELM_FULL_KERNEL(15)
DEFINE_AX_HELM_FULL_KERNEL(16)


__kernel
void ax_helm_stress_kernel_vector_part2(__global real * __restrict__ au,
                                        __global real * __restrict__ av,
                                        __global real * __restrict__ aw,
                                        __global const real * __restrict__ u,
                                        __global const real * __restrict__ v,
                                        __global const real * __restrict__ w,
                                        __global const real * __restrict__ h2,
                                        __global const real * __restrict__ B,
                                        const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  for (int i = idx; i < n; i += str) {
    au[i] = au[i] + h2[i] * B[i] * u[i];
    av[i] = av[i] + h2[i] * B[i] * v[i];
    aw[i] = aw[i] + h2[i] * B[i] * w[i];
  }
  
}

#endif // __MATH_AX_HELM_FULL_KERNEL_CL__
