#ifndef __MATH_AX_HELM_SVV_FACTORIZED_KERNEL_H__
#define __MATH_AX_HELM_SVV_FACTORIZED_KERNEL_H__
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
 * Fused device kernel for the factorized SVV Helmholtz operator
 *
 * Evaluates, in a single pass over each element,
 *   w = D^T (h1 G) D u + B^T (h1_svv G) B u,
 * with B = (Br, Bs, Bt) the complementary (filtered) derivative matrices.
 * This mirrors ax_helm_svv_factorized_lx in the CPU backend: both
 * contractions share the k-plane sweep, the geometric factors are loaded
 * once, and w is stored once (no read-modify-write of the output).
 *
 * The kernel is kept separate from the autotuned ax_helm kernels on purpose,
 * so that those keep their register footprint untouched.
 *
 * PAD selects a shared memory pitch of LX+1, removing bank conflicts when LX
 * is a power of two.
 */
template< typename T, const int LX, const bool PAD >
__global__ void __launch_bounds__(LX*LX)
ax_helm_svv_factorized_kernel(T * __restrict__ w,
                              const T * __restrict__ u,
                              const T * __restrict__ dx,
                              const T * __restrict__ dy,
                              const T * __restrict__ dz,
                              const T * __restrict__ h1,
                              const T * __restrict__ br,
                              const T * __restrict__ bs,
                              const T * __restrict__ bt,
                              const T * __restrict__ h1_svv,
                              const T * __restrict__ g11,
                              const T * __restrict__ g22,
                              const T * __restrict__ g33,
                              const T * __restrict__ g12,
                              const T * __restrict__ g13,
                              const T * __restrict__ g23) {

  constexpr int P = PAD ? LX + 1 : LX;

  __shared__ T shdx[LX * P];
  __shared__ T shdy[LX * P];
  __shared__ T shdz[LX * P];
  __shared__ T shbr[LX * P];
  __shared__ T shbs[LX * P];
  __shared__ T shbt[LX * P];

  __shared__ T shu[LX * P];
  __shared__ T shur[LX * LX];      // only accessed along fastest dimension
  __shared__ T shus[LX * P];
  __shared__ T shur_svv[LX * LX];  // only accessed along fastest dimension
  __shared__ T shus_svv[LX * P];

  T ru[LX];
  T rw[LX];
  T rut;
  T rut_svv;

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ij_p = i + j*P;
  const int ele = e*LX*LX*LX;

  shdx[ij_p] = dx[ij];
  shdy[ij_p] = dy[ij];
  shdz[ij_p] = dz[ij];
  shbr[ij_p] = br[ij];
  shbs[ij_p] = bs[ij];
  shbt[ij_p] = bt[ij];

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k*LX*LX + ele];
    rw[k] = 0.0;
  }

  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    const T G00 = g11[ijk+ele];
    const T G11 = g22[ijk+ele];
    const T G22 = g33[ijk+ele];
    const T G01 = g12[ijk+ele];
    const T G02 = g13[ijk+ele];
    const T G12 = g23[ijk+ele];

    T ttmp = 0.0;
    T ttmp_svv = 0.0;
    shu[ij_p] = ru[k];
#pragma unroll
    for (int l = 0; l < LX; l++) {
      ttmp += shdz[k+l*P] * ru[l];
      ttmp_svv += shbt[k+l*P] * ru[l];
    }
    __syncthreads();

    T rtmp = 0.0;
    T stmp = 0.0;
    T rtmp_svv = 0.0;
    T stmp_svv = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      const T ul_r = shu[l+j*P];
      const T ul_s = shu[i+l*P];
      rtmp += shdx[i+l*P] * ul_r;
      stmp += shdy[j+l*P] * ul_s;
      rtmp_svv += shbr[i+l*P] * ul_r;
      stmp_svv += shbs[j+l*P] * ul_s;
    }

    const T H1 = h1[ijk+ele];
    shur[ij] = H1 * (G00 * rtmp + G01 * stmp + G02 * ttmp);
    shus[ij_p] = H1 * (G01 * rtmp + G11 * stmp + G12 * ttmp);
    rut = H1 * (G02 * rtmp + G12 * stmp + G22 * ttmp);

    const T H1_SVV = h1_svv[ijk+ele];
    shur_svv[ij] = H1_SVV
                 * (G00 * rtmp_svv + G01 * stmp_svv + G02 * ttmp_svv);
    shus_svv[ij_p] = H1_SVV
                   * (G01 * rtmp_svv + G11 * stmp_svv + G12 * ttmp_svv);
    rut_svv = H1_SVV
            * (G02 * rtmp_svv + G12 * stmp_svv + G22 * ttmp_svv);

    __syncthreads();

    T wijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      wijke += shur[l+j*LX] * shdx[l+i*P]
             + shur_svv[l+j*LX] * shbr[l+i*P];
      wijke += shus[i+l*P] * shdy[l+j*P]
             + shus_svv[i+l*P] * shbs[l+j*P];
      rw[l] += rut * shdz[k+l*P] + rut_svv * shbt[k+l*P];
    }
    rw[k] += wijke;
  }

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    w[ij + k*LX*LX + ele] = rw[k];
  }
}

#endif // __MATH_AX_HELM_SVV_FACTORIZED_KERNEL_H__
