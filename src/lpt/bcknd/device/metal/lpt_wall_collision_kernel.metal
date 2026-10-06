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
 * Metal compute kernels for elastic particle-wall collisions.
 *
 * Metal allows at most 31 buffer arguments per function, which the
 * single CUDA/HIP kernel exceeds. The work is therefore split: the lag
 * and laglag vector triples are reflected by their own kernel, and the
 * main kernel reflects the current state and writes the new position.
 * Hit detection depends only on the pre-update position, so the lag
 * kernels must run before the main kernel, which overwrites it.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Facet geometry helpers. Indices follow the Fortran convention and are
 * 1-based; the facet normals are stored as (a, b, facet, element).
 */

static int lpt_wall_dof_idx(const int i, const int j, const int k,
                            const int e, const int lx, const int ly,
                            const int lz) {
  return (i - 1) + lx * ((j - 1) + ly * ((k - 1) + lz * (e - 1)));
}

static int lpt_wall_mask_idx(const int facet, const int el,
                             const int n_facets) {
  return (facet - 1) + n_facets * (el - 1);
}

static int lpt_wall_imax(const int a, const int b) {
  return a > b ? a : b;
}

static float lpt_norm3(const float x, const float y, const float z) {
  return sqrt(x * x + y * y + z * z);
}

/** Facet normal at a query point, see coef_get_normal */
static void lpt_coef_get_normal(device const float *nx,
                                device const float *ny,
                                device const float *nz,
                                const int i, const int j, const int k,
                                const int e, const int facet, const int lx,
                                thread float *normal_x,
                                thread float *normal_y,
                                thread float *normal_z) {
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
    *normal_x = 0.0f;
    *normal_y = 0.0f;
    *normal_z = 0.0f;
    return;
  }

  if (a < 1 || a > lx || b < 1 || b > lx || e < 1) {
    *normal_x = 0.0f;
    *normal_y = 0.0f;
    *normal_z = 0.0f;
    return;
  }

  const int normal_idx = (a - 1) + lx * ((b - 1) + lx *
                                         ((facet - 1) + 6 * (e - 1)));

  *normal_x = nx[normal_idx];
  *normal_y = ny[normal_idx];
  *normal_z = nz[normal_idx];
}

static void lpt_wall_facet_center(device const float *dm_x,
                                  device const float *dm_y,
                                  device const float *dm_z,
                                  const int lx, const int ly, const int lz,
                                  const int el, const int facet,
                                  thread float *wall_x,
                                  thread float *wall_y,
                                  thread float *wall_z) {
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
    *wall_x = 0.0f;
    *wall_y = 0.0f;
    *wall_z = 0.0f;
    return;
  }

  *wall_x = dm_x[idx];
  *wall_y = dm_y[idx];
  *wall_z = dm_z[idx];
}

static void lpt_wall_facet_normal(device const float *nx,
                                  device const float *ny,
                                  device const float *nz,
                                  const int lx, const int ly, const int lz,
                                  const int el, const int facet,
                                  thread float *normal_x,
                                  thread float *normal_y,
                                  thread float *normal_z) {
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
    *normal_x = 0.0f;
    *normal_y = 0.0f;
    *normal_z = 0.0f;
    break;
  }
}

static float lpt_signed_plane_distance(const float x, const float y,
                                       const float z, const float wall_x,
                                       const float wall_y, const float wall_z,
                                       const float normal_x,
                                       const float normal_y,
                                       const float normal_z) {
  const float eps = 1.0e-12f;
  const float nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return -1.0e30f;
  }

  return ((x - wall_x) * normal_x + (y - wall_y) * normal_y +
          (z - wall_z) * normal_z) / nmag;
}

static bool lpt_wall_facet_is_hit(device const int *wall_facet_mask,
                                  device const float *dm_x,
                                  device const float *dm_y,
                                  device const float *dm_z,
                                  device const float *nx,
                                  device const float *ny,
                                  device const float *nz,
                                  const int lx, const int ly, const int lz,
                                  const int el, const int facet,
                                  const int gdim,
                                  const float x_old, const float y_old,
                                  const float z_old,
                                  const float x, const float y,
                                  const float z, const float radius) {
  const int n_facets = 2 * gdim;
  const float tol = 1.0e-8f;
  float normal_x, normal_y, normal_z;
  float wall_x, wall_y, wall_z;

  if (facet < 1 || facet > n_facets) {
    return false;
  }
  if (wall_facet_mask[lpt_wall_mask_idx(facet, el, n_facets)] == 0) {
    return false;
  }

  lpt_wall_facet_normal(nx, ny, nz, lx, ly, lz, el, facet, &normal_x,
                        &normal_y, &normal_z);
  if (lpt_norm3(normal_x, normal_y, normal_z) <= 1.0e-12f) {
    return false;
  }

  lpt_wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el, facet, &wall_x,
                        &wall_y, &wall_z);
  const float dist_old = lpt_signed_plane_distance(x_old, y_old, z_old,
                                                   wall_x, wall_y, wall_z,
                                                   normal_x, normal_y,
                                                   normal_z);
  const float dist_new = lpt_signed_plane_distance(x, y, z, wall_x, wall_y,
                                                   wall_z, normal_x,
                                                   normal_y, normal_z);
  const float penetration = dist_new + radius;

  if (penetration <= tol) {
    return false;
  }
  if (dist_new <= dist_old + tol) {
    return false;
  }

  return true;
}

static void lpt_reflect_position(thread float *x, thread float *y,
                                 thread float *z,
                                 const float wall_x, const float wall_y,
                                 const float wall_z, const float normal_x,
                                 const float normal_y, const float normal_z,
                                 const float radius) {
  const float eps = 1.0e-12f;
  const float nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return;
  }

  const float nhat_x = normal_x / nmag;
  const float nhat_y = normal_y / nmag;
  const float nhat_z = normal_z / nmag;
  const float signed_contact_distance =
    (*x - wall_x) * nhat_x + (*y - wall_y) * nhat_y +
    (*z - wall_z) * nhat_z + radius;

  if (signed_contact_distance <= 0.0f) {
    return;
  }

  *x -= 2.0f * signed_contact_distance * nhat_x;
  *y -= 2.0f * signed_contact_distance * nhat_y;
  *z -= 2.0f * signed_contact_distance * nhat_z;
}

/** Reflect one vector triple stored at entry @a p of three arrays */
static void lpt_reflect_triple(device float *ax, device float *ay,
                               device float *az, const int p,
                               const float normal_x, const float normal_y,
                               const float normal_z) {
  const float eps = 1.0e-12f;
  const float nmag = lpt_norm3(normal_x, normal_y, normal_z);
  if (nmag <= eps) {
    return;
  }

  const float nhat_x = normal_x / nmag;
  const float nhat_y = normal_y / nmag;
  const float nhat_z = normal_z / nmag;
  const float vn = ax[p] * nhat_x + ay[p] * nhat_y + az[p] * nhat_z;

  ax[p] -= 2.0f * vn * nhat_x;
  ay[p] -= 2.0f * vn * nhat_y;
  az[p] -= 2.0f * vn * nhat_z;
}

/**
 * Reflect the current velocity, previous velocity and acceleration of
 * every colliding particle, and write the reflected position.
 */
kernel void lpt_handle_elastic_wall_collisions_kernel(
    device const int *wall_facet_mask [[ buffer(0) ]],
    device const int *el_list [[ buffer(1) ]],
    device const float *x_old [[ buffer(2) ]],
    device const float *y_old [[ buffer(3) ]],
    device const float *z_old [[ buffer(4) ]],
    device float *x [[ buffer(5) ]],
    device float *y [[ buffer(6) ]],
    device float *z [[ buffer(7) ]],
    device const float *d [[ buffer(8) ]],
    device float *u [[ buffer(9) ]],
    device float *v [[ buffer(10) ]],
    device float *w [[ buffer(11) ]],
    device float *u_old [[ buffer(12) ]],
    device float *v_old [[ buffer(13) ]],
    device float *w_old [[ buffer(14) ]],
    device float *acc_x [[ buffer(15) ]],
    device float *acc_y [[ buffer(16) ]],
    device float *acc_z [[ buffer(17) ]],
    device const float *dm_x [[ buffer(18) ]],
    device const float *dm_y [[ buffer(19) ]],
    device const float *dm_z [[ buffer(20) ]],
    device const float *nx [[ buffer(21) ]],
    device const float *ny [[ buffer(22) ]],
    device const float *nz [[ buffer(23) ]],
    constant int &n [[ buffer(24) ]],
    constant int &gdim [[ buffer(25) ]],
    constant int &nelv [[ buffer(26) ]],
    constant int &lx [[ buffer(27) ]],
    constant int &ly [[ buffer(28) ]],
    constant int &lz [[ buffer(29) ]],
    uint gid [[ thread_position_in_grid ]]) {
  if (gid >= (uint) n) return;

  const int p = (int) gid;
  const int n_facets = 2 * gdim;

  const int el = el_list[p];
  if (el < 0) {
    return;
  }

  const int el_mesh = el + 1;
  if (el_mesh > nelv) {
    return;
  }

  const float radius = 0.5f * d[p];
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
    return;
  }

  float x_new = x[p];
  float y_new = y[p];
  float z_new = z[p];

  for (int hit_idx = 0; hit_idx < hit_count; hit_idx++) {
    const int facet = hit_facets[hit_idx];
    float normal_x, normal_y, normal_z;
    float wall_x, wall_y, wall_z;

    lpt_wall_facet_normal(nx, ny, nz, lx, ly, lz, el_mesh, facet,
                          &normal_x, &normal_y, &normal_z);
    if (lpt_norm3(normal_x, normal_y, normal_z) <= 1.0e-12f) {
      continue;
    }

    lpt_wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el_mesh, facet,
                          &wall_x, &wall_y, &wall_z);
    lpt_reflect_position(&x_new, &y_new, &z_new, wall_x, wall_y, wall_z,
                         normal_x, normal_y, normal_z, radius);
    lpt_reflect_triple(u, v, w, p, normal_x, normal_y, normal_z);
    lpt_reflect_triple(u_old, v_old, w_old, p, normal_x, normal_y, normal_z);
    lpt_reflect_triple(acc_x, acc_y, acc_z, p, normal_x, normal_y, normal_z);
  }

  x[p] = x_new;
  y[p] = y_new;
  z[p] = z_new;
}

/**
 * Reflect one lagged velocity triple and one lagged acceleration triple.
 * Detection repeats the main kernel's, so this must be dispatched before
 * the main kernel overwrites the position.
 */
kernel void lpt_wall_collision_reflect_lag_kernel(
    device const int *wall_facet_mask [[ buffer(0) ]],
    device const int *el_list [[ buffer(1) ]],
    device const float *x_old [[ buffer(2) ]],
    device const float *y_old [[ buffer(3) ]],
    device const float *z_old [[ buffer(4) ]],
    device const float *x [[ buffer(5) ]],
    device const float *y [[ buffer(6) ]],
    device const float *z [[ buffer(7) ]],
    device const float *d [[ buffer(8) ]],
    device float *u_lag [[ buffer(9) ]],
    device float *v_lag [[ buffer(10) ]],
    device float *w_lag [[ buffer(11) ]],
    device float *acc_xlag [[ buffer(12) ]],
    device float *acc_ylag [[ buffer(13) ]],
    device float *acc_zlag [[ buffer(14) ]],
    device const float *dm_x [[ buffer(15) ]],
    device const float *dm_y [[ buffer(16) ]],
    device const float *dm_z [[ buffer(17) ]],
    device const float *nx [[ buffer(18) ]],
    device const float *ny [[ buffer(19) ]],
    device const float *nz [[ buffer(20) ]],
    constant int &n [[ buffer(21) ]],
    constant int &gdim [[ buffer(22) ]],
    constant int &nelv [[ buffer(23) ]],
    constant int &lx [[ buffer(24) ]],
    constant int &ly [[ buffer(25) ]],
    constant int &lz [[ buffer(26) ]],
    uint gid [[ thread_position_in_grid ]]) {
  if (gid >= (uint) n) return;

  const int p = (int) gid;
  const int n_facets = 2 * gdim;

  const int el = el_list[p];
  if (el < 0) {
    return;
  }

  const int el_mesh = el + 1;
  if (el_mesh > nelv) {
    return;
  }

  const float radius = 0.5f * d[p];

  for (int candidate = 1; candidate <= n_facets; candidate++) {
    if (!lpt_wall_facet_is_hit(wall_facet_mask, dm_x, dm_y, dm_z, nx, ny,
                               nz, lx, ly, lz, el_mesh, candidate, gdim,
                               x_old[p], y_old[p], z_old[p], x[p], y[p],
                               z[p], radius)) {
      continue;
    }

    float normal_x, normal_y, normal_z;
    lpt_wall_facet_normal(nx, ny, nz, lx, ly, lz, el_mesh, candidate,
                          &normal_x, &normal_y, &normal_z);
    if (lpt_norm3(normal_x, normal_y, normal_z) <= 1.0e-12f) {
      continue;
    }

    lpt_reflect_triple(u_lag, v_lag, w_lag, p, normal_x, normal_y, normal_z);
    lpt_reflect_triple(acc_xlag, acc_ylag, acc_zlag, p, normal_x, normal_y,
                       normal_z);
  }
}
