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

#ifndef __LPT_WALL_COLLISION_KERNEL_H__
#define __LPT_WALL_COLLISION_KERNEL_H__

#include <math.h>
#include <sem/bcknd/device/hip/coef_kernel.h>

__device__ inline int lpt_wall_dof_idx(const int i,
                                       const int j,
                                       const int k,
                                       const int e,
                                       const int lx,
                                       const int ly,
                                       const int lz) {
  return (i - 1) + lx * ((j - 1) + ly * ((k - 1) + lz * (e - 1)));
}

__device__ inline int lpt_wall_mask_idx(const int facet,
                                        const int el,
                                        const int n_facets) {
  return (facet - 1) + n_facets * (el - 1);
}

__device__ inline int lpt_wall_imax(const int a, const int b) {
  return a > b ? a : b;
}

template<typename T>
__device__ inline T lpt_norm3(const T x, const T y, const T z) {
  return sqrt(x * x + y * y + z * z);
}

template<typename T>
__device__ inline void lpt_wall_facet_center(
    const T * __restrict__ dm_x,
    const T * __restrict__ dm_y,
    const T * __restrict__ dm_z,
    const int lx,
    const int ly,
    const int lz,
    const int el,
    const int facet,
    T &wall_x,
    T &wall_y,
    T &wall_z) {
  const int ic = lpt_wall_imax(1, (lx + 1) / 2);
  const int jc = lpt_wall_imax(1, (ly + 1) / 2);
  const int kc = lpt_wall_imax(1, (lz + 1) / 2);
  int idx = 0;

  switch (facet) {
  case 1:
    idx = lpt_wall_dof_idx(1, jc, kc, el, lx, ly, lz);
    break;
  case 2:
    idx = lpt_wall_dof_idx(lx, jc, kc, el, lx, ly, lz);
    break;
  case 3:
    idx = lpt_wall_dof_idx(ic, 1, kc, el, lx, ly, lz);
    break;
  case 4:
    idx = lpt_wall_dof_idx(ic, ly, kc, el, lx, ly, lz);
    break;
  case 5:
    idx = lpt_wall_dof_idx(ic, jc, 1, el, lx, ly, lz);
    break;
  case 6:
    idx = lpt_wall_dof_idx(ic, jc, lz, el, lx, ly, lz);
    break;
  default:
    wall_x = 0.0;
    wall_y = 0.0;
    wall_z = 0.0;
    return;
  }

  wall_x = dm_x[idx];
  wall_y = dm_y[idx];
  wall_z = dm_z[idx];
}

template<typename T>
__device__ inline void lpt_wall_facet_normal(
    const T * __restrict__ nx,
    const T * __restrict__ ny,
    const T * __restrict__ nz,
    const int lx,
    const int ly,
    const int lz,
    const int el,
    const int facet,
    T &normal_x,
    T &normal_y,
    T &normal_z) {
  const int ic = lpt_wall_imax(1, (lx + 1) / 2);
  const int jc = lpt_wall_imax(1, (ly + 1) / 2);
  const int kc = lpt_wall_imax(1, (lz + 1) / 2);

  switch (facet) {
  case 1:
  case 2:
    coef_get_normal_device(nx, ny, nz, 1, jc, kc, el, facet, lx, normal_x,
                           normal_y, normal_z);
    break;
  case 3:
  case 4:
    coef_get_normal_device(nx, ny, nz, ic, 1, kc, el, facet, lx, normal_x,
                           normal_y, normal_z);
    break;
  case 5:
  case 6:
    coef_get_normal_device(nx, ny, nz, ic, jc, 1, el, facet, lx, normal_x,
                           normal_y, normal_z);
    break;
  default:
    normal_x = 0.0;
    normal_y = 0.0;
    normal_z = 0.0;
    break;
  }
}

template<typename T>
__device__ inline T lpt_signed_plane_distance(const T x,
                                              const T y,
                                              const T z,
                                              const T wall_x,
                                              const T wall_y,
                                              const T wall_z,
                                              const T normal_x,
                                              const T normal_y,
                                              const T normal_z) {
  const T eps = (T) 1.0e-12;
  const T nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return -1.0e30;
  }

  return ((x - wall_x) * normal_x + (y - wall_y) * normal_y +
          (z - wall_z) * normal_z) / nmag;
}

template<typename T>
__device__ inline bool lpt_wall_facet_is_hit(
    const int * __restrict__ wall_facet_mask,
    const T * __restrict__ dm_x,
    const T * __restrict__ dm_y,
    const T * __restrict__ dm_z,
    const T * __restrict__ nx,
    const T * __restrict__ ny,
    const T * __restrict__ nz,
    const int lx,
    const int ly,
    const int lz,
    const int el,
    const int facet,
    const int gdim,
    const T x_old,
    const T y_old,
    const T z_old,
    const T x,
    const T y,
    const T z,
    const T radius) {
  const int n_facets = 2 * gdim;
  const T tol = (T) 1.0e-8;
  T normal_x, normal_y, normal_z;
  T wall_x, wall_y, wall_z;

  if (facet < 1 || facet > n_facets) {
    return false;
  }
  if (wall_facet_mask[lpt_wall_mask_idx(facet, el, n_facets)] == 0) {
    return false;
  }

  lpt_wall_facet_normal(nx, ny, nz, lx, ly, lz, el, facet, normal_x,
                        normal_y, normal_z);
  if (lpt_norm3(normal_x, normal_y, normal_z) <= (T) 1.0e-12) {
    return false;
  }

  lpt_wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el, facet, wall_x,
                        wall_y, wall_z);
  const T dist_old = lpt_signed_plane_distance(x_old, y_old, z_old, wall_x,
                                               wall_y, wall_z, normal_x,
                                               normal_y, normal_z);
  const T dist_new = lpt_signed_plane_distance(x, y, z, wall_x, wall_y,
                                               wall_z, normal_x, normal_y,
                                               normal_z);
  const T penetration = dist_new + radius;

  if (penetration <= tol) {
    return false;
  }
  if (dist_new <= dist_old + tol) {
    return false;
  }

  return true;
}

template<typename T>
__device__ inline void lpt_reflect_position(T &x,
                                            T &y,
                                            T &z,
                                            const T wall_x,
                                            const T wall_y,
                                            const T wall_z,
                                            const T normal_x,
                                            const T normal_y,
                                            const T normal_z,
                                            const T radius) {
  const T eps = (T) 1.0e-12;
  const T nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return;
  }

  const T nhat_x = normal_x / nmag;
  const T nhat_y = normal_y / nmag;
  const T nhat_z = normal_z / nmag;
  const T signed_contact_distance =
      (x - wall_x) * nhat_x + (y - wall_y) * nhat_y +
      (z - wall_z) * nhat_z + radius;

  if (signed_contact_distance <= 0.0) {
    return;
  }

  x -= 2.0 * signed_contact_distance * nhat_x;
  y -= 2.0 * signed_contact_distance * nhat_y;
  z -= 2.0 * signed_contact_distance * nhat_z;
}

template<typename T>
__device__ inline void lpt_reflect_vector_components(T &x,
                                                     T &y,
                                                     T &z,
                                                     const T normal_x,
                                                     const T normal_y,
                                                     const T normal_z) {
  const T eps = (T) 1.0e-12;
  const T nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return;
  }

  const T nhat_x = normal_x / nmag;
  const T nhat_y = normal_y / nmag;
  const T nhat_z = normal_z / nmag;
  const T vn = x * nhat_x + y * nhat_y + z * nhat_z;

  x -= 2.0 * vn * nhat_x;
  y -= 2.0 * vn * nhat_y;
  z -= 2.0 * vn * nhat_z;
}

template<typename T>
__global__ void lpt_handle_elastic_wall_collisions_kernel(
    const int * __restrict__ wall_facet_mask,
    const int * __restrict__ el_list,
    const T * __restrict__ x_old,
    const T * __restrict__ y_old,
    const T * __restrict__ z_old,
    T * __restrict__ x,
    T * __restrict__ y,
    T * __restrict__ z,
    const T * __restrict__ d,
    T * __restrict__ u,
    T * __restrict__ v,
    T * __restrict__ w,
    T * __restrict__ u_lag,
    T * __restrict__ v_lag,
    T * __restrict__ w_lag,
    T * __restrict__ u_laglag,
    T * __restrict__ v_laglag,
    T * __restrict__ w_laglag,
    T * __restrict__ acc_xlag,
    T * __restrict__ acc_ylag,
    T * __restrict__ acc_zlag,
    T * __restrict__ acc_xlaglag,
    T * __restrict__ acc_ylaglag,
    T * __restrict__ acc_zlaglag,
    T * __restrict__ u_old,
    T * __restrict__ v_old,
    T * __restrict__ w_old,
    T * __restrict__ acc_x,
    T * __restrict__ acc_y,
    T * __restrict__ acc_z,
    const T * __restrict__ dm_x,
    const T * __restrict__ dm_y,
    const T * __restrict__ dm_z,
    const T * __restrict__ nx,
    const T * __restrict__ ny,
    const T * __restrict__ nz,
    const int n,
    const int gdim,
    const int nelv,
    const int lx,
    const int ly,
    const int lz,
    const int lag_len) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const int n_facets = 2 * gdim;

  for (int p = idx; p < n; p += str) {
    const int el = el_list[p];
    if (el < 0) {
      continue;
    }

    const int el_mesh = el + 1;
    if (el_mesh > nelv) {
      continue;
    }

    const T radius = 0.5 * d[p];
    int hit_facets[6];
    int hit_count = 0;

    for (int candidate = 1; candidate <= n_facets; candidate++) {
      if (lpt_wall_facet_is_hit(wall_facet_mask, dm_x, dm_y, dm_z, nx, ny,
                                nz, lx, ly, lz, el_mesh, candidate, gdim,
                                x_old[p], y_old[p], z_old[p], x[p], y[p],
                                z[p], radius)) {
        hit_facets[hit_count] = candidate;
        hit_count++;
      }
    }

    if (hit_count == 0) {
      continue;
    }

    T x_new = x[p];
    T y_new = y[p];
    T z_new = z[p];

    for (int hit_idx = 0; hit_idx < hit_count; hit_idx++) {
      const int facet = hit_facets[hit_idx];
      T normal_x, normal_y, normal_z;
      T wall_x, wall_y, wall_z;

      lpt_wall_facet_normal(nx, ny, nz, lx, ly, lz, el_mesh, facet,
                            normal_x, normal_y, normal_z);
      if (lpt_norm3(normal_x, normal_y, normal_z) <= (T) 1.0e-12) {
        continue;
      }

      lpt_wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el_mesh, facet,
                            wall_x, wall_y, wall_z);
      lpt_reflect_position(x_new, y_new, z_new, wall_x, wall_y, wall_z,
                           normal_x, normal_y, normal_z, radius);
      lpt_reflect_vector_components(u[p], v[p], w[p], normal_x, normal_y,
                                    normal_z);
      lpt_reflect_vector_components(u_old[p], v_old[p], w_old[p], normal_x,
                                    normal_y, normal_z);
      lpt_reflect_vector_components(acc_x[p], acc_y[p], acc_z[p], normal_x,
                                    normal_y, normal_z);

      if (lag_len >= 1) {
        lpt_reflect_vector_components(u_lag[p], v_lag[p], w_lag[p],
                                      normal_x, normal_y, normal_z);
        lpt_reflect_vector_components(acc_xlag[p], acc_ylag[p], acc_zlag[p],
                                      normal_x, normal_y, normal_z);
      }

      if (lag_len >= 2) {
        lpt_reflect_vector_components(u_laglag[p], v_laglag[p], w_laglag[p],
                                      normal_x, normal_y, normal_z);
        lpt_reflect_vector_components(acc_xlaglag[p], acc_ylaglag[p],
                                      acc_zlaglag[p], normal_x, normal_y,
                                      normal_z);
      }
    }

    x[p] = x_new;
    y[p] = y_new;
    z[p] = z_new;
  }
}

#endif
