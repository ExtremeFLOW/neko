#ifndef __SEM_COEF_KERNEL_H__
#define __SEM_COEF_KERNEL_H__
/*
 Copyright (c) 2022-2026, The Neko Authors
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
 * Device kernel for coef geometry
 */
template< typename T, const int LX, const int CHUNKS >
__global__ void coef_generate_geo_kernel(T * __restrict__ G11,
					 T * __restrict__ G12,
					 T * __restrict__ G13,
					 T * __restrict__ G22,
					 T * __restrict__ G23,
					 T * __restrict__ G33,
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
					 const int gdim) {

  int i,j,k;

  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  __shared__ T shw3[LX * LX * LX];

  j = iii;
  while( j < (LX * LX * LX)) {
    const int i = j + e * LX * LX * LX;
    G11[i] = (drdx[i]*drdx[i] + drdy[i]*drdy[i] + drdz[i]*drdz[i]) * jacinv[i];
    G22[i] = (dsdx[i]*dsdx[i] + dsdy[i]*dsdy[i] + dsdz[i]*dsdz[i]) * jacinv[i];
    G33[i] = (dtdx[i]*dtdx[i] + dtdy[i]*dtdy[i] + dtdz[i]*dtdz[i]) * jacinv[i];

    G12[i] = (drdx[i]*dsdx[i] + drdy[i]*dsdy[i] + drdz[i]*dsdz[i]) * jacinv[i];
    G13[i] = (drdx[i]*dtdx[i] + drdy[i]*dtdy[i] + drdz[i]*dtdz[i]) * jacinv[i];
    G23[i] = (dsdx[i]*dtdx[i] + dsdy[i]*dtdy[i] + dsdz[i]*dtdz[i]) * jacinv[i];

    shw3[j] = w3[j];

    j = j + CHUNKS;
  }

  __syncthreads();

  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      G11[ijk + e * LX * LX * LX] *= shw3[ijk];
      G12[ijk + e * LX * LX * LX] *= shw3[ijk];
      G13[ijk + e * LX * LX * LX] *= shw3[ijk];
      G22[ijk + e * LX * LX * LX] *= shw3[ijk];
      G23[ijk + e * LX * LX * LX] *= shw3[ijk];
      G33[ijk + e * LX * LX * LX] *= shw3[ijk];
    }
  }
}

/**
 * Device kernel for coef dxyz
 */
template< typename T, const int LX, const int CHUNKS >
__global__ void coef_generate_dxyz_kernel(T * __restrict__ dxdr,
					  T * __restrict__ dydr,
					  T * __restrict__ dzdr,
					  T * __restrict__ dxds,
					  T * __restrict__ dyds,
					  T * __restrict__ dzds,
					  T * __restrict__ dxdt,
					  T * __restrict__ dydt,
					  T * __restrict__ dzdt,
					  const T * __restrict__ dx,
					  const T * __restrict__ dy,
					  const T * __restrict__ dz,
					  const T * __restrict__ x,
					  const T * __restrict__ y,
					  const T * __restrict__ z) {

  int i,j,k;

  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  __shared__ T shu[LX * LX * LX];

  if (iii < (LX * LX)) {
    shdx[iii] = dx[iii];
    shdy[iii] = dy[iii];
    shdz[iii] = dz[iii];
  }

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = x[j + e * LX * LX * LX];
    j = j + CHUNKS;
  }

  __syncthreads();

  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];
      }
      dxdr[ijk + e * LX * LX * LX] = rtmp;
      dxds[ijk + e * LX * LX * LX] = stmp;
      dxdt[ijk + e * LX * LX * LX] = ttmp;
    }
  }

  __syncthreads();

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = y[j + e * LX * LX * LX];
    j = j + CHUNKS;
  }

  __syncthreads();

  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];
      }
      dydr[ijk + e * LX * LX * LX] = rtmp;
      dyds[ijk + e * LX * LX * LX] = stmp;
      dydt[ijk + e * LX * LX * LX] = ttmp;
    }
  }

  __syncthreads();

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = z[j + e * LX * LX * LX];
    j = j + CHUNKS;
  }

  __syncthreads();

  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];
      }
      dzdr[ijk + e * LX * LX * LX] = rtmp;
      dzds[ijk + e * LX * LX * LX] = stmp;
      dzdt[ijk + e * LX * LX * LX] = ttmp;
    }
  }
}

/**
 * Device kernel for coef drst
 */
template< typename T >
__global__ void coef_generate_drst_kernel(T * __restrict__ jac,
					  T * __restrict__ jacinv,
					  T * __restrict__ drdx,
					  T * __restrict__ drdy,
					  T * __restrict__ drdz,
					  T * __restrict__ dsdx,
					  T * __restrict__ dsdy,
					  T * __restrict__ dsdz,
					  T * __restrict__ dtdx,
					  T * __restrict__ dtdy,
					  T * __restrict__ dtdz,
					  const T * __restrict__ dxdr,
					  const T * __restrict__ dydr,
					  const T * __restrict__ dzdr,
					  const T * __restrict__ dxds,
					  const T * __restrict__ dyds,
					  const T * __restrict__ dzds,
					  const T * __restrict__ dxdt,
					  const T * __restrict__ dydt,
					  const T * __restrict__ dzdt,
					  const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const T one = 1.0;

  for (int i = idx; i < n; i += str) {
    jac[i] = (dxdr[i] * dyds[i] * dzdt[i])
           + (dxdt[i] * dydr[i] * dzds[i])
           + (dxds[i] * dydt[i] * dzdr[i])
           - (dxdr[i] * dydt[i] * dzds[i])
           - (dxds[i] * dydr[i] * dzdt[i])
           - (dxdt[i] * dyds[i] * dzdr[i]);
    jacinv[i] = one / jac[i];

    drdx[i] = dyds[i]*dzdt[i] - dydt[i]*dzds[i];
    drdy[i] = dxdt[i]*dzds[i] - dxds[i]*dzdt[i];
    drdz[i] = dxds[i]*dydt[i] - dxdt[i]*dyds[i];
    dsdx[i] = dydt[i]*dzdr[i] - dydr[i]*dzdt[i];
    dsdy[i] = dxdr[i]*dzdt[i] - dxdt[i]*dzdr[i];
    dsdz[i] = dxdt[i]*dydr[i] - dxdr[i]*dydt[i];
    dtdx[i] = dydr[i]*dzds[i] - dyds[i]*dzdr[i];
    dtdy[i] = dxds[i]*dzdr[i] - dxdr[i]*dzds[i];
    dtdz[i] = dxdr[i]*dyds[i] - dxds[i]*dydr[i];

  }

}

/**
 * Device kernel for coef_generate_mass
 */
template <typename T>
__global__ void coef_generate_mass_kernel(T * __restrict__ B,
                                          T * __restrict__ Binv,
                                          const T * __restrict__ jac,
                                          const T * __restrict__ w3,
                                          int lxyz, int nel) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int n = lxyz * nel;

  if (idx < n) {
    int local_idx = idx - (idx / lxyz) * lxyz;

    T mass_val = jac[idx] * w3[local_idx];

    B[idx] = mass_val;
    Binv[idx] = mass_val;
  }
}

/**
 * Device kernel for coef_generate_area_and_normal
 */
template< typename T, const int LX >
__global__
void coef_generate_area_and_normal_kernel(T * __restrict__ area,
                                          T * __restrict__ nx,
                                          T * __restrict__ ny,
                                          T * __restrict__ nz,
                                          const T * __restrict__ dxdr,
                                          const T * __restrict__ dydr,
                                          const T * __restrict__ dzdr,
                                          const T * __restrict__ dxds,
                                          const T * __restrict__ dyds,
                                          const T * __restrict__ dzds,
                                          const T * __restrict__ dxdt,
                                          const T * __restrict__ dydt,
                                          const T * __restrict__ dzdt,
                                          const T * __restrict__ wx,
                                          const T * __restrict__ wy,
                                          const T * __restrict__ wz,
                                          const T eps) {
  int i, j, k;
  int f, out_idx;
  const T one = 1.0;
  const T m_one = -1.0;
  T tx, ty, tz, dot, weight, length, sgn;

  const int e = blockIdx.x;
  const int iii = threadIdx.x;

  const int lxyz = LX * LX * LX;
  const int lxy = LX * LX;

  for (int ijk = iii; ijk < lxyz; ijk += blockDim.x) {

    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;

    const int offset = ijk + (e * lxyz);
    const int face_offset = e * lxy * 6;

    // ds x dt
    if (i == 0 || i == LX - 1) {
      tx = dyds[offset] * dzdt[offset] - dzds[offset] * dydt[offset];
      ty = dzds[offset] * dxdt[offset] - dxds[offset] * dzdt[offset];
      tz = dxds[offset] * dydt[offset] - dyds[offset] * dxdt[offset];

      dot = tx*tx + ty*ty + tz*tz;
      length = sqrt(dot);
      weight = wy[j] * wz[k];

      if (i == 0) {
        f = 0;
        sgn = m_one;
      } else {
        f = 1;
        sgn = one;
      }

      out_idx = j + (k * LX) + (f * lxy) + face_offset;

      area[out_idx] = length * weight;

      if (length > eps) {
        nx[out_idx] = (tx / length) * sgn;
        ny[out_idx] = (ty / length) * sgn;
        nz[out_idx] = (tz / length) * sgn;
      } else {
        nx[out_idx] = tx * sgn;
        ny[out_idx] = ty * sgn;
        nz[out_idx] = tz * sgn;
      }
    }

    // dr x dt
    if (j == 0 || j == LX - 1) {
      tx = dydr[offset] * dzdt[offset] - dzdr[offset] * dydt[offset];
      ty = dzdr[offset] * dxdt[offset] - dxdr[offset] * dzdt[offset];
      tz = dxdr[offset] * dydt[offset] - dydr[offset] * dxdt[offset];

      dot = tx*tx + ty*ty + tz*tz;
      length = sqrt(dot);
      weight = wx[i] * wz[k];

      if (j == 0) {
        f = 2;
        sgn = one;
      } else {
        f = 3;
        sgn = m_one;
      }

      out_idx = i + (k * LX) + (f * lxy) + face_offset;

      area[out_idx] = length * weight;

      if (length > eps) {
        nx[out_idx] = (tx / length) * sgn;
        ny[out_idx] = (ty / length) * sgn;
        nz[out_idx] = (tz / length) * sgn;
      } else {
        nx[out_idx] = tx * sgn;
        ny[out_idx] = ty * sgn;
        nz[out_idx] = tz * sgn;
      }
    }

    // dr x ds
    if (k == 0 || k == LX - 1) {
      tx = dydr[offset] * dzds[offset] - dzdr[offset] * dyds[offset];
      ty = dzdr[offset] * dxds[offset] - dxdr[offset] * dzds[offset];
      tz = dxdr[offset] * dyds[offset] - dydr[offset] * dxds[offset];

      dot = tx*tx + ty*ty + tz*tz;
      length = sqrt(dot);
      weight = wx[i] * wy[j];

      if (k == 0) {
        f = 4;
        sgn = m_one;
      } else {
        f = 5;
        sgn = one;
      }

      out_idx = i + (j * LX) + (f * lxy) + face_offset;

      area[out_idx] = length * weight;

      if (length > eps) {
        nx[out_idx] = (tx / length) * sgn;
        ny[out_idx] = (ty / length) * sgn;
        nz[out_idx] = (tz / length) * sgn;
      } else {
        nx[out_idx] = tx * sgn;
        ny[out_idx] = ty * sgn;
        nz[out_idx] = tz * sgn;
      }
    }
  }
}


template< typename T >
__device__ inline void coef_get_normal_device(const T * __restrict__ nx,
                                              const T * __restrict__ ny,
                                              const T * __restrict__ nz,
                                              const int i,
                                              const int j,
                                              const int k,
                                              const int e,
                                              const int facet,
                                              const int lx,
                                              T &normal_x,
                                              T &normal_y,
                                              T &normal_z) {
  int a = 0;
  int b = 0;

  switch (facet) {
  case 1:
  case 2:
    a = j;
    b = k;
    break;
  case 3:
  case 4:
    a = i;
    b = k;
    break;
  case 5:
  case 6:
    a = i;
    b = j;
    break;
  default:
    normal_x = 0.0;
    normal_y = 0.0;
    normal_z = 0.0;
    return;
  }

  if (a < 1 || a > lx || b < 1 || b > lx || e < 1) {
    normal_x = 0.0;
    normal_y = 0.0;
    normal_z = 0.0;
    return;
  }

  const int normal_idx = (a - 1) + lx * ((b - 1) + lx *
                         ((facet - 1) + 6 * (e - 1)));

  normal_x = nx[normal_idx];
  normal_y = ny[normal_idx];
  normal_z = nz[normal_idx];
}

/**
 * Device kernel for coef_get_normal.
 *
 * Query indices use the Fortran coef_get_normal convention: i, j, k, e and
 * facet are 1-based. Output arrays are 0-based device vectors of length n.
 */
template< typename T >
__global__
void coef_get_normal_kernel(T * __restrict__ normal_x,
                            T * __restrict__ normal_y,
                            T * __restrict__ normal_z,
                            const T * __restrict__ nx,
                            const T * __restrict__ ny,
                            const T * __restrict__ nz,
                            const int * __restrict__ i_idx,
                            const int * __restrict__ j_idx,
                            const int * __restrict__ k_idx,
                            const int * __restrict__ e_idx,
                            const int * __restrict__ facet_idx,
                            const int lx,
                            const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int p = idx; p < n; p += str) {
    coef_get_normal_device(nx, ny, nz, i_idx[p], j_idx[p], k_idx[p], e_idx[p],
                           facet_idx[p], lx, normal_x[p], normal_y[p],
                           normal_z[p]);
  }
}

#endif // __SEM_COEF_KERNEL_H__
