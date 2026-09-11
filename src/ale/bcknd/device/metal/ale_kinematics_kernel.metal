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
 * Metal compute kernels for ALE mesh kinematics.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Rigid body kinematics, passed by value.
 * @note The layout must match kinematics_params_t in the CUDA, HIP and
 * OpenCL backends and the bind(c) type in ale_routines_device.F90.
 */
struct kinematics_params_t {
  float cx, cy, cz;
  float vtx, vty, vtz;
  float vax, vay, vaz;
  float px, py, pz;
  float r11, r12, r13;
  float r21, r22, r23;
  float r31, r32, r33;
};

kernel void ale_add_kinematics_kernel(constant int &n [[ buffer(0) ]],
                                      device float *wx [[ buffer(1) ]],
                                      device float *wy [[ buffer(2) ]],
                                      device float *wz [[ buffer(3) ]],
                                      device const float *x_ref [[ buffer(4) ]],
                                      device const float *y_ref [[ buffer(5) ]],
                                      device const float *z_ref [[ buffer(6) ]],
                                      device const float *phi [[ buffer(7) ]],
                                      device const float *x [[ buffer(8) ]],
                                      device const float *y [[ buffer(9) ]],
                                      device const float *z [[ buffer(10) ]],
                                      constant kinematics_params_t &kin_params [[ buffer(11) ]],
                                      uint gid [[ thread_position_in_grid ]]) {
  if (gid >= (uint) n) return;

  const int i = (int) gid;
  const float tol = 1e-6f;

  const float p_val = phi[i];
  float v_tan_x, v_tan_y, v_tan_z;

  if (fabs(p_val - 1.0f) < tol) {
    const float rx = x[i] - kin_params.cx;
    const float ry = y[i] - kin_params.cy;
    const float rz = z[i] - kin_params.cz;

    v_tan_x = kin_params.vay * rz - kin_params.vaz * ry;
    v_tan_y = kin_params.vaz * rx - kin_params.vax * rz;
    v_tan_z = kin_params.vax * ry - kin_params.vay * rx;
  } else {
    const float dx_ref = x_ref[i] - kin_params.px;
    const float dy_ref = y_ref[i] - kin_params.py;
    const float dz_ref = z_ref[i] - kin_params.pz;

    const float rx_target = kin_params.r11 * dx_ref +
      kin_params.r12 * dy_ref +
      kin_params.r13 * dz_ref;

    const float ry_target = kin_params.r21 * dx_ref +
      kin_params.r22 * dy_ref +
      kin_params.r23 * dz_ref;

    const float rz_target = kin_params.r31 * dx_ref +
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

kernel void compute_cheap_dist_kernel(device float *d [[ buffer(0) ]],
                                      device const float *x [[ buffer(1) ]],
                                      device const float *y [[ buffer(2) ]],
                                      device const float *z [[ buffer(3) ]],
                                      constant int &lx [[ buffer(4) ]],
                                      constant int &ly [[ buffer(5) ]],
                                      constant int &lz [[ buffer(6) ]],
                                      constant int &nel [[ buffer(7) ]],
                                      constant int &local_iters [[ buffer(8) ]],
                                      device atomic_int *nchange [[ buffer(9) ]],
                                      uint gid [[ thread_position_in_grid ]]) {
  if (gid >= (uint) nel) return;

  const int e = (int) gid;

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
          const float x1 = x[idx1];
          const float y1 = y[idx1];
          const float z1 = z[idx1];
          float d1 = d[idx1];

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
                const float x2 = x[idx2];
                const float y2 = y[idx2];
                const float z2 = z[idx2];
                const float d2 = d[idx2];

                const float dist = sqrt((x1 - x2) * (x1 - x2) +
                                        (y1 - y2) * (y1 - y2) +
                                        (z1 - z2) * (z1 - z2));
                const float dtmp = d2 + dist;

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
    atomic_fetch_add_explicit(nchange, 1, memory_order_relaxed);
  }
}
