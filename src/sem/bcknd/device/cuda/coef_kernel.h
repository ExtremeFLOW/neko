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

/*
 * Elements per thread block for the coef kstep kernels
 *
 * A kstep kernel gives one (LX,LX) thread plane to an element and walks the
 * element in k, so a block holds LX*LX*EB threads for the EB elements stacked
 * along threadIdx.z. This replaces the earlier layout of 1024 threads per
 * element, under which the product of block size and per thread register count
 * could exceed the 65536 registers a block can be given. Built for sm_120 with
 * nvcc 13 in double precision, the area and normal kernel takes 70 registers
 * per thread at LX = 2 and the dxyz kernel 126 at LX = 16, so at 1024 threads
 * both failed to launch with "too many resources requested for launch".
 *
 * EB is fixed at compile time rather than tuned. The SEM operator kernels pick
 * theirs with an autotuner (see math/bcknd/device/cuda/elem_block.h), but the
 * coef kernels run once per mesh, and once per step only when the mesh moves,
 * so there is nothing for a tuning sweep to amortise against. Taking the
 * largest power of two that keeps a block at or below
 * COEF_EB_MAX_BLOCK_THREADS threads puts every LX in 2..16 between 64 and 256
 * threads per block.
 */
#define COEF_EB_MAX_BLOCK_THREADS 256

template< int LX >
struct coef_elem_block {
  static const int raw = COEF_EB_MAX_BLOCK_THREADS / (LX * LX);
  static const int value = (raw >= 16) ? 16 : (raw >= 8) ? 8 :
                           (raw >= 4) ? 4 : (raw >= 2) ? 2 : 1;
};

#define COEF_EB(LX) (coef_elem_block<LX>::value)
#define COEF_EB_NTHRDS(LX) dim3((LX), (LX), COEF_EB(LX))
#define COEF_EB_NBLCKS(NELV, LX)                                               \
  dim3(((NELV) + COEF_EB(LX) - 1) / COEF_EB(LX), 1, 1)

/**
 * Device kernel for coef geometry
 *
 * Each thread owns one (i,j) column of an element and walks it in k. The
 * quadrature weights w3 are the same for every element, so they are read from
 * global memory at the point where they are applied; the earlier version
 * staged them in a LX*LX*LX shared array and wrote G11..G33 twice, once
 * without the weight and once with it.
 */
template< typename T, const int LX, const int EB >
__global__ void __launch_bounds__(LX * LX * EB)
coef_generate_geo_kernel(T * __restrict__ G11,
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
                         const int nelv,
                         const int gdim) {

  const int e = blockIdx.x * EB + ((EB == 1) ? 0 : threadIdx.z);

  /* The kernel has no barriers, so the threads of a partially filled tail
     block can simply leave */
  if (e >= nelv) {
    return;
  }

  const int ij = threadIdx.x + threadIdx.y * LX;
  const int ele = e * LX * LX * LX;

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k * LX * LX;
    const int idx = ijk + ele;

    const T w = w3[ijk];
    const T jinv = jacinv[idx];

    const T rx = drdx[idx];
    const T ry = drdy[idx];
    const T rz = drdz[idx];
    const T sx = dsdx[idx];
    const T sy = dsdy[idx];
    const T sz = dsdz[idx];
    const T tx = dtdx[idx];
    const T ty = dtdy[idx];
    const T tz = dtdz[idx];

    G11[idx] = (rx*rx + ry*ry + rz*rz) * jinv * w;
    G22[idx] = (sx*sx + sy*sy + sz*sz) * jinv * w;
    G33[idx] = (tx*tx + ty*ty + tz*tz) * jinv * w;

    G12[idx] = (rx*sx + ry*sy + rz*sz) * jinv * w;
    G13[idx] = (rx*tx + ry*ty + rz*tz) * jinv * w;
    G23[idx] = (sx*tx + sy*ty + sz*tz) * jinv * w;
  }
}

/**
 * Compute the r, s and t derivatives of one coordinate of an element
 *
 * Called once per coordinate by coef_generate_dxyz_kernel with the thread
 * layout of that kernel: the caller's thread owns the (i,j) column given by
 * ij, sh is the offset of the calling element's plane in shu, and ele is the
 * offset of the element in the coordinate and derivative arrays. Threads of a
 * partially filled tail block have active false; they still take part in the
 * barriers and only their stores are dropped.
 */
template< typename T, const int LX >
__device__ inline void coef_dxyz_component(T * __restrict__ dudr,
                                           T * __restrict__ duds,
                                           T * __restrict__ dudt,
                                           const T * __restrict__ u,
                                           T * shu,
                                           const T * shdx,
                                           const T * shdy,
                                           const T * shdz,
                                           const int i,
                                           const int j,
                                           const int sh,
                                           const int ele,
                                           const bool active) {
  const int ij = i + j * LX;
  T ru[LX];

#pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k * LX * LX + ele];
  }

  __syncthreads();

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k * LX * LX;

    /* The t derivative runs along the column this thread holds in ru */
    T ttmp = 0.0;
    shu[sh + ij] = ru[k];
#pragma unroll
    for (int l = 0; l < LX; l++) {
      ttmp += shdz[k + l * LX] * ru[l];
    }

    __syncthreads();

    /* The r and s derivatives run along the k plane now in shu */
    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      rtmp += shdx[i + l * LX] * shu[sh + l + j * LX];
      stmp += shdy[j + l * LX] * shu[sh + i + l * LX];
    }

    if (active) {
      dudr[ijk + ele] = rtmp;
      duds[ijk + ele] = stmp;
      dudt[ijk + ele] = ttmp;
    }

    __syncthreads();
  }
}

/**
 * Device kernel for coef dxyz
 *
 * Each thread owns one (i,j) column of an element and walks it in k, holding
 * the column of the coordinate being differentiated in registers for the t
 * derivative and staging one k plane of it in shared memory for the r and s
 * derivatives. The three coordinates are taken in turn through the same
 * shared plane.
 */
template< typename T, const int LX, const int EB >
__global__ void __launch_bounds__(LX * LX * EB)
coef_generate_dxyz_kernel(T * __restrict__ dxdr,
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
                          const T * __restrict__ z,
                          const int nelv) {

  /* Element independent, one copy per block */
  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  /* One k plane per element in the block */
  __shared__ T shu[EB * LX * LX];

  const int eb = (EB == 1) ? 0 : threadIdx.z;
  const int e_blk = blockIdx.x * EB + eb;
  /* Threads past the last element still have to reach the barriers in the k
     loop, so clamp their reads and drop their stores rather than returning
     early. At EB == 1 this all constant folds away */
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int sh = eb * LX * LX;
  const int i = threadIdx.x;
  const int j = threadIdx.y;
  const int ij = i + j * LX;
  const int ele = e * LX * LX * LX;

  if (eb == 0) {
    shdx[ij] = dx[ij];
    shdy[ij] = dy[ij];
    shdz[ij] = dz[ij];
  }

  coef_dxyz_component<T, LX>(dxdr, dxds, dxdt, x,
                             shu, shdx, shdy, shdz, i, j, sh, ele, active);
  coef_dxyz_component<T, LX>(dydr, dyds, dydt, y,
                             shu, shdx, shdy, shdz, i, j, sh, ele, active);
  coef_dxyz_component<T, LX>(dzdr, dzds, dzdt, z,
                             shu, shdx, shdy, shdz, i, j, sh, ele, active);
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
 *
 * Each thread owns one (i,j) column of an element and walks it in k, writing
 * the facet entries its column lands on: the r facets 0 and 1 when i is 0 or
 * LX-1, the s facets 2 and 3 when j is, and the t facets 4 and 5 when k is.
 */
template< typename T, const int LX, const int EB >
__global__ void __launch_bounds__(LX * LX * EB)
coef_generate_area_and_normal_kernel(T * __restrict__ area,
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
                                     const T eps,
                                     const int nelv) {
  int f, out_idx;
  const T one = 1.0;
  const T m_one = -1.0;
  T tx, ty, tz, dot, weight, length, sgn;

  const int e = blockIdx.x * EB + ((EB == 1) ? 0 : threadIdx.z);

  /* The kernel has no barriers, so the threads of a partially filled tail
     block can simply leave */
  if (e >= nelv) {
    return;
  }

  const int i = threadIdx.x;
  const int j = threadIdx.y;

  const int lxyz = LX * LX * LX;
  const int lxy = LX * LX;

  const int face_offset = e * lxy * 6;
  const T wxi = wx[i];
  const T wyj = wy[j];

  for (int k = 0; k < LX; ++k) {

    const int ijk = i + (j * LX) + (k * lxy);
    const int offset = ijk + (e * lxyz);

    // ds x dt
    if (i == 0 || i == LX - 1) {
      tx = dyds[offset] * dzdt[offset] - dzds[offset] * dydt[offset];
      ty = dzds[offset] * dxdt[offset] - dxds[offset] * dzdt[offset];
      tz = dxds[offset] * dydt[offset] - dyds[offset] * dxdt[offset];

      dot = tx*tx + ty*ty + tz*tz;
      length = sqrt(dot);
      weight = wyj * wz[k];

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
      weight = wxi * wz[k];

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
      weight = wxi * wyj;

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
