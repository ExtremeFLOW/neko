#ifndef __ALE_KINEMATICS_KERNEL_CL__
#define __ALE_KINEMATICS_KERNEL_CL__
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
 * Rigid body kinematics, passed by value.
 * @note The layout must match kinematics_params_t in the CUDA, HIP and
 * Metal backends and the bind(c) type in ale_routines_device.F90.
 */
typedef struct {
  real cx, cy, cz;
  real vtx, vty, vtz;
  real vax, vay, vaz;
  real px, py, pz;
  real r11, r12, r13;
  real r21, r22, r23;
  real r31, r32, r33;
} kinematics_params_t;

/** Device kernel for add_kinematics_to_mesh_velocity */
__kernel void ale_add_kinematics_kernel(const int n,
                                        __global real * __restrict__ wx,
                                        __global real * __restrict__ wy,
                                        __global real * __restrict__ wz,
                                        __global const real * __restrict__ x_ref,
                                        __global const real * __restrict__ y_ref,
                                        __global const real * __restrict__ z_ref,
                                        __global const real * __restrict__ phi,
                                        __global const real * __restrict__ x,
                                        __global const real * __restrict__ y,
                                        __global const real * __restrict__ z,
                                        const kinematics_params_t kin_params) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  const real one = (real) 1.0;
  const real tol = (real) 1e-6;

  for (int i = idx; i < n; i += str) {
    const real p_val = phi[i];
    real v_tan_x, v_tan_y, v_tan_z;

    if (fabs(p_val - one) < tol) {
      const real rx = x[i] - kin_params.cx;
      const real ry = y[i] - kin_params.cy;
      const real rz = z[i] - kin_params.cz;

      v_tan_x = kin_params.vay * rz - kin_params.vaz * ry;
      v_tan_y = kin_params.vaz * rx - kin_params.vax * rz;
      v_tan_z = kin_params.vax * ry - kin_params.vay * rx;
    } else {
      const real dx_ref = x_ref[i] - kin_params.px;
      const real dy_ref = y_ref[i] - kin_params.py;
      const real dz_ref = z_ref[i] - kin_params.pz;

      const real rx_target = kin_params.r11 * dx_ref +
        kin_params.r12 * dy_ref +
        kin_params.r13 * dz_ref;

      const real ry_target = kin_params.r21 * dx_ref +
        kin_params.r22 * dy_ref +
        kin_params.r23 * dz_ref;

      const real rz_target = kin_params.r31 * dx_ref +
        kin_params.r32 * dy_ref +
        kin_params.r33 * dz_ref;

      v_tan_x = kin_params.vay * rz_target - kin_params.vaz * ry_target;
      v_tan_y = kin_params.vaz * rx_target - kin_params.vax * rz_target;
      v_tan_z = kin_params.vax * ry_target - kin_params.vay * rx_target;
    }

    wx[i] += (kin_params.vtx + v_tan_x) * p_val;
    wy[i] += (kin_params.vty + v_tan_y) * p_val;
    wz[i] += (kin_params.vtz + v_tan_z) * p_val;
  }
}

/** Device kernel for compute_cheap_dist */
__kernel void compute_cheap_dist_kernel(__global real * __restrict__ d,
                                        __global const real * __restrict__ x,
                                        __global const real * __restrict__ y,
                                        __global const real * __restrict__ z,
                                        const int lx, const int ly,
                                        const int lz, const int nel,
                                        const int local_iters,
                                        __global int * __restrict__ nchange) {

  const int e = get_global_id(0);
  if (e >= nel) return;

  int iter = 1;
  bool element_changed_ever = false;
  bool changed_local = true;

  const int lxy = lx * ly;
  const int lxyz = lx * ly * lz;
  const int e_offset = e * lxyz;

  while (changed_local && iter <= local_iters) {
    changed_local = false;

    /* Loop over GLL nodes in this element */
    for (int k = 0; k < lz; ++k) {
      for (int j = 0; j < ly; ++j) {
        for (int i = 0; i < lx; ++i) {
          const int idx1 = i + j * lx + k * lxy + e_offset;
          const real x1 = x[idx1];
          const real y1 = y[idx1];
          const real z1 = z[idx1];
          real d1 = d[idx1];

          const int i0 = max(0, i - 1);
          const int i1 = min(lx - 1, i + 1);
          const int j0 = max(0, j - 1);
          const int j1 = min(ly - 1, j + 1);
          const int k0 = max(0, k - 1);
          const int k1 = min(lz - 1, k + 1);

          /* Neighbor check */
          for (int kk = k0; kk <= k1; ++kk) {
            for (int jj = j0; jj <= j1; ++jj) {
              for (int ii = i0; ii <= i1; ++ii) {
                if (ii == i && jj == j && kk == k) continue;

                const int idx2 = ii + jj * lx + kk * lxy + e_offset;
                const real x2 = x[idx2];
                const real y2 = y[idx2];
                const real z2 = z[idx2];
                const real d2 = d[idx2];

                const real dist = sqrt((x1 - x2) * (x1 - x2) +
                                       (y1 - y2) * (y1 - y2) +
                                       (z1 - z2) * (z1 - z2));
                const real dtmp = d2 + dist;

                if (dtmp < d1) {
                  d1 = dtmp;
                  d[idx1] = d1;   /* Update locally */
                  changed_local = true;
                }
              }
            }
          }
        }
      }
    }
    if (changed_local) element_changed_ever = true;
    iter++;
  }

  if (element_changed_ever) {
    atomic_add(nchange, 1);
  }
}

#endif // __ALE_KINEMATICS_KERNEL_CL__
