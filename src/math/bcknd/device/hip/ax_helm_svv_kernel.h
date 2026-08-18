#ifndef __MATH_AX_HELM_SVV_KERNEL_H__
#define __MATH_AX_HELM_SVV_KERNEL_H__
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
 * Fused device kernel for the asymmetric SVV Helmholtz operator.
 *
 * The filter matrices are selected by the caller. An identity matrix in a
 * direction disables filtering in that direction while retaining the complete
 * Helmholtz gradient and divergence.
 */
template<typename T, const int LX>
__global__ void ax_helm_svv_kernel(
    T * __restrict__ w,
    const T * __restrict__ u,
    const T * __restrict__ dx,
    const T * __restrict__ dy,
    const T * __restrict__ dz,
    const T * __restrict__ h1,
    const T * __restrict__ drdx,
    const T * __restrict__ drdy,
    const T * __restrict__ drdz,
    const T * __restrict__ dsdx,
    const T * __restrict__ dsdy,
    const T * __restrict__ dsdz,
    const T * __restrict__ dtdx,
    const T * __restrict__ dtdy,
    const T * __restrict__ dtdz,
    const T * __restrict__ jacinv,
    const T * __restrict__ w3,
    const T * __restrict__ h1_svv,
    const T * __restrict__ filter_r,
    const T * __restrict__ filter_s,
    const T * __restrict__ filter_t) {

  extern __shared__ T shared[];
  T *shfield = shared;
  T *shwork = shared + LX * LX * LX;

  const int e = blockIdx.x;
  const int i = threadIdx.x;
  const int j = threadIdx.y;
  const int ij = i + j * LX;
  const int lx2 = LX * LX;
  const int elem = e * LX * lx2;

#pragma unroll 1
  for (int k = 0; k < LX; ++k) {
    w[ij + k * lx2 + elem] = 0.0;
  }

  // Process the x, y and z components of the physical gradient in turn.
#pragma unroll 1
  for (int component = 0; component < 3; ++component) {

    // Form one component of the physical gradient.
#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T ur = 0.0;
      T us = 0.0;
      T ut = 0.0;

#pragma unroll
      for (int l = 0; l < LX; ++l) {
        ur += dx[i + l * LX] * u[l + j * LX + k * lx2 + elem];
        us += dy[j + l * LX] * u[i + l * LX + k * lx2 + elem];
        ut += dz[k + l * LX] * u[ij + l * lx2 + elem];
      }

      const int ijk = ij + k * lx2;
      const int index = ijk + elem;
      if (component == 0) {
        shwork[ijk] = (ur * drdx[index] + us * dsdx[index] +
                       ut * dtdx[index]) * jacinv[index];
      }
      else if (component == 1) {
        shwork[ijk] = (ur * drdy[index] + us * dsdy[index] +
                       ut * dtdy[index]) * jacinv[index];
      }
      else {
        shwork[ijk] = (ur * drdz[index] + us * dsdz[index] +
                       ut * dtdz[index]) * jacinv[index];
      }
    }

    // Apply filter_r, filter_s and filter_t as a tensor product.
    __syncthreads();
#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T value = 0.0;
#pragma unroll
      for (int l = 0; l < LX; ++l) {
        value += filter_r[i + l * LX] *
                 shwork[l + j * LX + k * lx2];
      }
      shfield[ij + k * lx2] = value;
    }
    __syncthreads();

#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T value = 0.0;
#pragma unroll
      for (int l = 0; l < LX; ++l) {
        value += filter_s[l + j * LX] *
                 shfield[i + l * LX + k * lx2];
      }
      shwork[ij + k * lx2] = value;
    }
    __syncthreads();

#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T value = 0.0;
#pragma unroll
      for (int l = 0; l < LX; ++l) {
        value += filter_t[l + k * LX] * shwork[ij + l * lx2];
      }
      shfield[ij + k * lx2] = value;
    }
    __syncthreads();

    // Recompute the unfiltered gradient and form the weighted physical flux.
#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T ur = 0.0;
      T us = 0.0;
      T ut = 0.0;
#pragma unroll
      for (int l = 0; l < LX; ++l) {
        ur += dx[i + l * LX] * u[l + j * LX + k * lx2 + elem];
        us += dy[j + l * LX] * u[i + l * LX + k * lx2 + elem];
        ut += dz[k + l * LX] * u[ij + l * lx2 + elem];
      }

      const int ijk = ij + k * lx2;
      const int index = ijk + elem;
      T gradient;
      if (component == 0) {
        gradient = (ur * drdx[index] + us * dsdx[index] +
                    ut * dtdx[index]) * jacinv[index];
      }
      else if (component == 1) {
        gradient = (ur * drdy[index] + us * dsdy[index] +
                    ut * dtdy[index]) * jacinv[index];
      }
      else {
        gradient = (ur * drdz[index] + us * dsdz[index] +
                    ut * dtdz[index]) * jacinv[index];
      }

      const T weighted_gradient =
          w3[ijk] * (h1[index] * gradient +
                     h1_svv[index] * (gradient - shfield[ijk]));
      shwork[ijk] = weighted_gradient;
    }
    __syncthreads();

    // Apply the r-direction reference flux contribution.
#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      const int ijk = ij + k * lx2;
      const int index = ijk + elem;
      if (component == 0) {
        shfield[ijk] = drdx[index] * shwork[ijk];
      }
      else if (component == 1) {
        shfield[ijk] = drdy[index] * shwork[ijk];
      }
      else {
        shfield[ijk] = drdz[index] * shwork[ijk];
      }
    }
    __syncthreads();

#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T value = 0.0;
#pragma unroll
      for (int l = 0; l < LX; ++l) {
        value += dx[l + i * LX] * shfield[l + j * LX + k * lx2];
      }
      w[ij + k * lx2 + elem] += value;
    }
    __syncthreads();

    // Apply the s-direction reference flux contribution.
#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      const int ijk = ij + k * lx2;
      const int index = ijk + elem;
      if (component == 0) {
        shfield[ijk] = dsdx[index] * shwork[ijk];
      }
      else if (component == 1) {
        shfield[ijk] = dsdy[index] * shwork[ijk];
      }
      else {
        shfield[ijk] = dsdz[index] * shwork[ijk];
      }
    }
    __syncthreads();

#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T value = 0.0;
#pragma unroll
      for (int l = 0; l < LX; ++l) {
        value += dy[l + j * LX] * shfield[i + l * LX + k * lx2];
      }
      w[ij + k * lx2 + elem] += value;
    }
    __syncthreads();

    // Apply the t-direction reference flux contribution.
#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      const int ijk = ij + k * lx2;
      const int index = ijk + elem;
      if (component == 0) {
        shfield[ijk] = dtdx[index] * shwork[ijk];
      }
      else if (component == 1) {
        shfield[ijk] = dtdy[index] * shwork[ijk];
      }
      else {
        shfield[ijk] = dtdz[index] * shwork[ijk];
      }
    }
    __syncthreads();

#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T value = 0.0;
#pragma unroll
      for (int l = 0; l < LX; ++l) {
        value += dz[l + k * LX] * shfield[ij + l * lx2];
      }
      w[ij + k * lx2 + elem] += value;
    }
    __syncthreads();
  }
}

#endif
