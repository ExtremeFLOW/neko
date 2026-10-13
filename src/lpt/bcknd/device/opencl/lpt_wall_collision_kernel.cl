#ifndef __LPT_WALL_COLLISION_KERNEL_CL__
#define __LPT_WALL_COLLISION_KERNEL_CL__
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

/*
 * Facet geometry helpers. Indices follow the Fortran convention and are
 * 1-based; the facet normals are stored as (a, b, facet, element).
 */

inline int lpt_wall_dof_idx(const int i, const int j, const int k,
                            const int e, const int lx, const int ly,
                            const int lz) {
  return (i - 1) + lx * ((j - 1) + ly * ((k - 1) + lz * (e - 1)));
}

inline int lpt_wall_mask_idx(const int facet, const int el,
                             const int n_facets) {
  return (facet - 1) + n_facets * (el - 1);
}

inline int lpt_wall_imax(const int a, const int b) {
  return a > b ? a : b;
}

inline real lpt_norm3(const real x, const real y, const real z) {
  return sqrt(x * x + y * y + z * z);
}

/** Facet normal at a query point, see coef_get_normal */
inline void lpt_coef_get_normal(__global const real * __restrict__ nx,
                                __global const real * __restrict__ ny,
                                __global const real * __restrict__ nz,
                                const int i, const int j, const int k,
                                const int e, const int facet, const int lx,
                                real *normal_x, real *normal_y,
                                real *normal_z) {
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
    *normal_x = (real) 0.0;
    *normal_y = (real) 0.0;
    *normal_z = (real) 0.0;
    return;
  }

  if (a < 1 || a > lx || b < 1 || b > lx || e < 1) {
    *normal_x = (real) 0.0;
    *normal_y = (real) 0.0;
    *normal_z = (real) 0.0;
    return;
  }

  const int normal_idx = (a - 1) + lx * ((b - 1) + lx *
                                         ((facet - 1) + 6 * (e - 1)));

  *normal_x = nx[normal_idx];
  *normal_y = ny[normal_idx];
  *normal_z = nz[normal_idx];
}

inline void lpt_wall_facet_center(__global const real * __restrict__ dm_x,
                                  __global const real * __restrict__ dm_y,
                                  __global const real * __restrict__ dm_z,
                                  const int lx, const int ly, const int lz,
                                  const int el, const int facet,
                                  real *wall_x, real *wall_y, real *wall_z) {
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
    *wall_x = (real) 0.0;
    *wall_y = (real) 0.0;
    *wall_z = (real) 0.0;
    return;
  }

  *wall_x = dm_x[idx];
  *wall_y = dm_y[idx];
  *wall_z = dm_z[idx];
}

inline void lpt_wall_facet_normal(__global const real * __restrict__ nx,
                                  __global const real * __restrict__ ny,
                                  __global const real * __restrict__ nz,
                                  const int lx, const int ly, const int lz,
                                  const int el, const int facet,
                                  real *normal_x, real *normal_y,
                                  real *normal_z) {
  const int ic = lpt_wall_imax(1, (lx + 1) / 2);
  const int jc = lpt_wall_imax(1, (ly + 1) / 2);
  const int kc = lpt_wall_imax(1, (lz + 1) / 2);

  switch (facet) {
  case 1:
  case 2:
    lpt_coef_get_normal(nx, ny, nz, 1, jc, kc, el, facet, lx, normal_x,
                        normal_y, normal_z);
    break;
  case 3:
  case 4:
    lpt_coef_get_normal(nx, ny, nz, ic, 1, kc, el, facet, lx, normal_x,
                        normal_y, normal_z);
    break;
  case 5:
  case 6:
    lpt_coef_get_normal(nx, ny, nz, ic, jc, 1, el, facet, lx, normal_x,
                        normal_y, normal_z);
    break;
  default:
    *normal_x = (real) 0.0;
    *normal_y = (real) 0.0;
    *normal_z = (real) 0.0;
    break;
  }
}

inline real lpt_signed_plane_distance(const real x, const real y,
                                      const real z, const real wall_x,
                                      const real wall_y, const real wall_z,
                                      const real normal_x,
                                      const real normal_y,
                                      const real normal_z) {
  const real eps = (real) 1.0e-12;
  const real nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return (real) -1.0e30;
  }

  return ((x - wall_x) * normal_x + (y - wall_y) * normal_y +
          (z - wall_z) * normal_z) / nmag;
}

inline int lpt_wall_facet_is_hit(
    __global const int * __restrict__ wall_facet_mask,
    __global const real * __restrict__ dm_x,
    __global const real * __restrict__ dm_y,
    __global const real * __restrict__ dm_z,
    __global const real * __restrict__ nx,
    __global const real * __restrict__ ny,
    __global const real * __restrict__ nz,
    const int lx, const int ly, const int lz,
    const int el, const int facet, const int gdim,
    const real x_old, const real y_old, const real z_old,
    const real x, const real y, const real z, const real radius) {
  const int n_facets = 2 * gdim;
  const real tol = (real) 1.0e-8;
  real normal_x, normal_y, normal_z;
  real wall_x, wall_y, wall_z;

  if (facet < 1 || facet > n_facets) {
    return 0;
  }
  if (wall_facet_mask[lpt_wall_mask_idx(facet, el, n_facets)] == 0) {
    return 0;
  }

  lpt_wall_facet_normal(nx, ny, nz, lx, ly, lz, el, facet, &normal_x,
                        &normal_y, &normal_z);
  if (lpt_norm3(normal_x, normal_y, normal_z) <= (real) 1.0e-12) {
    return 0;
  }

  lpt_wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el, facet, &wall_x,
                        &wall_y, &wall_z);
  const real dist_old = lpt_signed_plane_distance(x_old, y_old, z_old, wall_x,
                                                  wall_y, wall_z, normal_x,
                                                  normal_y, normal_z);
  const real dist_new = lpt_signed_plane_distance(x, y, z, wall_x, wall_y,
                                                  wall_z, normal_x, normal_y,
                                                  normal_z);
  const real penetration = dist_new + radius;

  if (penetration <= tol) {
    return 0;
  }
  if (dist_new <= dist_old + tol) {
    return 0;
  }

  return 1;
}

inline void lpt_reflect_position(real *x, real *y, real *z,
                                 const real wall_x, const real wall_y,
                                 const real wall_z, const real normal_x,
                                 const real normal_y, const real normal_z,
                                 const real radius) {
  const real eps = (real) 1.0e-12;
  const real nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return;
  }

  const real nhat_x = normal_x / nmag;
  const real nhat_y = normal_y / nmag;
  const real nhat_z = normal_z / nmag;
  const real signed_contact_distance =
    (*x - wall_x) * nhat_x + (*y - wall_y) * nhat_y +
    (*z - wall_z) * nhat_z + radius;

  if (signed_contact_distance <= (real) 0.0) {
    return;
  }

  *x -= (real) 2.0 * signed_contact_distance * nhat_x;
  *y -= (real) 2.0 * signed_contact_distance * nhat_y;
  *z -= (real) 2.0 * signed_contact_distance * nhat_z;
}

/** Reflect one vector triple stored at entry @a p of three arrays */
inline void lpt_reflect_triple(__global real * __restrict__ ax,
                               __global real * __restrict__ ay,
                               __global real * __restrict__ az,
                               const int p,
                               const real normal_x, const real normal_y,
                               const real normal_z) {
  const real eps = (real) 1.0e-12;
  const real nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return;
  }

  const real nhat_x = normal_x / nmag;
  const real nhat_y = normal_y / nmag;
  const real nhat_z = normal_z / nmag;
  const real vn = ax[p] * nhat_x + ay[p] * nhat_y + az[p] * nhat_z;

  ax[p] -= (real) 2.0 * vn * nhat_x;
  ay[p] -= (real) 2.0 * vn * nhat_y;
  az[p] -= (real) 2.0 * vn * nhat_z;
}

/** Device kernel for elastic particle-wall collisions */
__kernel void lpt_handle_elastic_wall_collisions_kernel(
    __global const int * __restrict__ wall_facet_mask,
    __global const int * __restrict__ el_list,
    __global const real * __restrict__ x_old,
    __global const real * __restrict__ y_old,
    __global const real * __restrict__ z_old,
    __global real * __restrict__ x,
    __global real * __restrict__ y,
    __global real * __restrict__ z,
    __global const real * __restrict__ d,
    __global real * __restrict__ u,
    __global real * __restrict__ v,
    __global real * __restrict__ w,
    __global real * __restrict__ u_lag,
    __global real * __restrict__ v_lag,
    __global real * __restrict__ w_lag,
    __global real * __restrict__ u_laglag,
    __global real * __restrict__ v_laglag,
    __global real * __restrict__ w_laglag,
    __global real * __restrict__ acc_xlag,
    __global real * __restrict__ acc_ylag,
    __global real * __restrict__ acc_zlag,
    __global real * __restrict__ acc_xlaglag,
    __global real * __restrict__ acc_ylaglag,
    __global real * __restrict__ acc_zlaglag,
    __global real * __restrict__ u_old,
    __global real * __restrict__ v_old,
    __global real * __restrict__ w_old,
    __global real * __restrict__ acc_x,
    __global real * __restrict__ acc_y,
    __global real * __restrict__ acc_z,
    __global const real * __restrict__ dm_x,
    __global const real * __restrict__ dm_y,
    __global const real * __restrict__ dm_z,
    __global const real * __restrict__ nx,
    __global const real * __restrict__ ny,
    __global const real * __restrict__ nz,
    const int n,
    const int gdim,
    const int nelv,
    const int lx,
    const int ly,
    const int lz,
    const int lag_len) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);
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

    const real radius = (real) 0.5 * d[p];
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

    real x_new = x[p];
    real y_new = y[p];
    real z_new = z[p];

    for (int hit_idx = 0; hit_idx < hit_count; hit_idx++) {
      const int facet = hit_facets[hit_idx];
      real normal_x, normal_y, normal_z;
      real wall_x, wall_y, wall_z;

      lpt_wall_facet_normal(nx, ny, nz, lx, ly, lz, el_mesh, facet,
                            &normal_x, &normal_y, &normal_z);
      if (lpt_norm3(normal_x, normal_y, normal_z) <= (real) 1.0e-12) {
        continue;
      }

      lpt_wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el_mesh, facet,
                            &wall_x, &wall_y, &wall_z);
      lpt_reflect_position(&x_new, &y_new, &z_new, wall_x, wall_y, wall_z,
                           normal_x, normal_y, normal_z, radius);
      lpt_reflect_triple(u, v, w, p, normal_x, normal_y, normal_z);
      lpt_reflect_triple(u_old, v_old, w_old, p, normal_x, normal_y,
                         normal_z);
      lpt_reflect_triple(acc_x, acc_y, acc_z, p, normal_x, normal_y,
                         normal_z);

      if (lag_len >= 1) {
        lpt_reflect_triple(u_lag, v_lag, w_lag, p, normal_x, normal_y,
                           normal_z);
        lpt_reflect_triple(acc_xlag, acc_ylag, acc_zlag, p, normal_x,
                           normal_y, normal_z);
      }

      if (lag_len >= 2) {
        lpt_reflect_triple(u_laglag, v_laglag, w_laglag, p, normal_x,
                           normal_y, normal_z);
        lpt_reflect_triple(acc_xlaglag, acc_ylaglag, acc_zlaglag, p,
                           normal_x, normal_y, normal_z);
      }
    }

    x[p] = x_new;
    y[p] = y_new;
    z[p] = z_new;
  }
}

#endif // __LPT_WALL_COLLISION_KERNEL_CL__
