#ifndef __MATH_AX_HELM_SVV_FULL_KERNEL_H__
#define __MATH_AX_HELM_SVV_FULL_KERNEL_H__
/*
 Copyright (c) 2025-2026, The Neko Authors
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
 * Evaluate one physical derivative of a field at a nodal point.
 */
template<typename T, const int LX>
__device__ __forceinline__ T ax_helm_svv_full_derivative(
    const T * __restrict__ field,
    const int component,
    const int i,
    const int j,
    const int k,
    const int elem,
    const T * __restrict__ dx,
    const T * __restrict__ dy,
    const T * __restrict__ dz,
    const T * __restrict__ drdx,
    const T * __restrict__ drdy,
    const T * __restrict__ drdz,
    const T * __restrict__ dsdx,
    const T * __restrict__ dsdy,
    const T * __restrict__ dsdz,
    const T * __restrict__ dtdx,
    const T * __restrict__ dtdy,
    const T * __restrict__ dtdz,
    const T * __restrict__ jacinv) {

  const int lx2 = LX * LX;
  const int ij = i + j * LX;
  const int ijk = ij + k * lx2;
  const int index = ijk + elem;
  T ur = 0.0;
  T us = 0.0;
  T ut = 0.0;

#pragma unroll
  for (int l = 0; l < LX; ++l) {
    ur += dx[i + l * LX] * field[l + j * LX + k * lx2 + elem];
    us += dy[j + l * LX] * field[i + l * LX + k * lx2 + elem];
    ut += dz[k + l * LX] * field[ij + l * lx2 + elem];
  }

  if (component == 0) {
    return (ur * drdx[index] + us * dsdx[index] +
            ut * dtdx[index]) * jacinv[index];
  }
  if (component == 1) {
    return (ur * drdy[index] + us * dsdy[index] +
            ut * dtdy[index]) * jacinv[index];
  }
  return (ur * drdz[index] + us * dsdz[index] +
          ut * dtdz[index]) * jacinv[index];
}

/**
 * Fused device kernel for the asymmetric full-stress SVV Helmholtz operator.
 *
 * Each physical strain component is formed, filtered and consumed entirely
 * within the kernel. Only two element-sized shared arrays are required.
 */
template<typename T, const int LX>
__global__ void ax_helm_svv_full_kernel(
    T * __restrict__ au,
    T * __restrict__ av,
    T * __restrict__ aw,
    const T * __restrict__ u,
    const T * __restrict__ v,
    const T * __restrict__ w,
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
    const int index = ij + k * lx2 + elem;
    au[index] = 0.0;
    av[index] = 0.0;
    aw[index] = 0.0;
  }

  // Process s11, s22, s33, s12, s13 and s23 in turn.
#pragma unroll 1
  for (int stress = 0; stress < 6; ++stress) {

    // Form one component of grad(u) + grad(u)^T.
#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T value;
      if (stress == 0) {
        value = 2.0 * ax_helm_svv_full_derivative<T, LX>(
            u, 0, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else if (stress == 1) {
        value = 2.0 * ax_helm_svv_full_derivative<T, LX>(
            v, 1, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else if (stress == 2) {
        value = 2.0 * ax_helm_svv_full_derivative<T, LX>(
            w, 2, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else if (stress == 3) {
        value = ax_helm_svv_full_derivative<T, LX>(
            u, 1, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv)
          + ax_helm_svv_full_derivative<T, LX>(
            v, 0, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else if (stress == 4) {
        value = ax_helm_svv_full_derivative<T, LX>(
            u, 2, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv)
          + ax_helm_svv_full_derivative<T, LX>(
            w, 0, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else {
        value = ax_helm_svv_full_derivative<T, LX>(
            v, 2, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv)
          + ax_helm_svv_full_derivative<T, LX>(
            w, 1, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      shwork[ij + k * lx2] = value;
    }

    // Apply the selected tensor-product low-pass filter.
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

    // Recompute the unfiltered strain and combine both viscosities.
#pragma unroll 1
    for (int k = 0; k < LX; ++k) {
      T value;
      if (stress == 0) {
        value = 2.0 * ax_helm_svv_full_derivative<T, LX>(
            u, 0, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else if (stress == 1) {
        value = 2.0 * ax_helm_svv_full_derivative<T, LX>(
            v, 1, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else if (stress == 2) {
        value = 2.0 * ax_helm_svv_full_derivative<T, LX>(
            w, 2, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else if (stress == 3) {
        value = ax_helm_svv_full_derivative<T, LX>(
            u, 1, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv)
          + ax_helm_svv_full_derivative<T, LX>(
            v, 0, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else if (stress == 4) {
        value = ax_helm_svv_full_derivative<T, LX>(
            u, 2, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv)
          + ax_helm_svv_full_derivative<T, LX>(
            w, 0, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }
      else {
        value = ax_helm_svv_full_derivative<T, LX>(
            v, 2, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv)
          + ax_helm_svv_full_derivative<T, LX>(
            w, 1, i, j, k, elem, dx, dy, dz,
            drdx, drdy, drdz, dsdx, dsdy, dsdz,
            dtdx, dtdy, dtdz, jacinv);
      }

      const int ijk = ij + k * lx2;
      const int index = ijk + elem;
      shwork[ijk] =
          w3[ijk] * (h1[index] * value +
                     h1_svv[index] * (value - shfield[ijk]));
    }
    __syncthreads();

    // A diagonal strain contributes once; a shear strain contributes twice.
    const int contributions = stress < 3 ? 1 : 2;
#pragma unroll 1
    for (int contribution = 0;
         contribution < contributions; ++contribution) {
      int output;
      int physical;
      if (stress < 3) {
        output = stress;
        physical = stress;
      }
      else if (stress == 3) {
        output = contribution;
        physical = 1 - contribution;
      }
      else if (stress == 4) {
        output = contribution == 0 ? 0 : 2;
        physical = contribution == 0 ? 2 : 0;
      }
      else {
        output = contribution == 0 ? 1 : 2;
        physical = contribution == 0 ? 2 : 1;
      }

      T *result = output == 0 ? au : (output == 1 ? av : aw);

      // r-direction reference flux and divergence.
#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        const int ijk = ij + k * lx2;
        const int index = ijk + elem;
        T metric = physical == 0 ? drdx[index]
                   : (physical == 1 ? drdy[index] : drdz[index]);
        shfield[ijk] = metric * shwork[ijk];
      }
      __syncthreads();

#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        T value = 0.0;
#pragma unroll
        for (int l = 0; l < LX; ++l) {
          value += dx[l + i * LX] * shfield[l + j * LX + k * lx2];
        }
        result[ij + k * lx2 + elem] += value;
      }
      __syncthreads();

      // s-direction reference flux and divergence.
#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        const int ijk = ij + k * lx2;
        const int index = ijk + elem;
        T metric = physical == 0 ? dsdx[index]
                   : (physical == 1 ? dsdy[index] : dsdz[index]);
        shfield[ijk] = metric * shwork[ijk];
      }
      __syncthreads();

#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        T value = 0.0;
#pragma unroll
        for (int l = 0; l < LX; ++l) {
          value += dy[l + j * LX] * shfield[i + l * LX + k * lx2];
        }
        result[ij + k * lx2 + elem] += value;
      }
      __syncthreads();

      // t-direction reference flux and divergence.
#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        const int ijk = ij + k * lx2;
        const int index = ijk + elem;
        T metric = physical == 0 ? dtdx[index]
                   : (physical == 1 ? dtdy[index] : dtdz[index]);
        shfield[ijk] = metric * shwork[ijk];
      }
      __syncthreads();

#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        T value = 0.0;
#pragma unroll
        for (int l = 0; l < LX; ++l) {
          value += dz[l + k * LX] * shfield[ij + l * lx2];
        }
        result[ij + k * lx2 + elem] += value;
      }
      __syncthreads();
    }
  }
}

#endif
