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

#ifndef __COMMON_ALE_KINEMATICS_KERNEL_H__
#define __COMMON_ALE_KINEMATICS_KERNEL_H__

#ifndef KINEMATICS_PARAMS_T_DEFINED
#define KINEMATICS_PARAMS_T_DEFINED
typedef struct {
    real cx, cy, cz;
    real vtx, vty, vtz;
    real vax, vay, vaz;
    real px, py, pz;
    real r11, r12, r13;
    real r21, r22, r23;
    real r31, r32, r33;
} kinematics_params_t;
#endif

template< typename T >
__global__ void ale_add_kinematics_kernel(const int n, 
                T * __restrict__ wx, 
                T * __restrict__ wy, 
                T * __restrict__ wz,
                const T * __restrict__ x_ref, 
                const T * __restrict__ y_ref, 
                const T * __restrict__ z_ref,
                const T * __restrict__ phi, 
                const T * __restrict__ x, 
                const T * __restrict__ y, 
                const T * __restrict__ z,
                const kinematics_params_t kin_params) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const T p_val = phi[i];
    T v_tan_x, v_tan_y, v_tan_z;

    if (abs(p_val - 1.0) < 1e-6) {
      const T rx = x[i] - kin_params.cx;
      const T ry = y[i] - kin_params.cy;
      const T rz = z[i] - kin_params.cz;

      v_tan_x = kin_params.vay * rz - kin_params.vaz * ry;
      v_tan_y = kin_params.vaz * rx - kin_params.vax * rz;
      v_tan_z = kin_params.vax * ry - kin_params.vay * rx;
    } 
    else {
      const T dx_ref = x_ref[i] - kin_params.px;
      const T dy_ref = y_ref[i] - kin_params.py;
      const T dz_ref = z_ref[i] - kin_params.pz;

      const T rx_target = kin_params.r11*dx_ref + 
                          kin_params.r12*dy_ref + 
                          kin_params.r13*dz_ref;

      const T ry_target = kin_params.r21*dx_ref + 
                          kin_params.r22*dy_ref + 
                          kin_params.r23*dz_ref;

      const T rz_target = kin_params.r31*dx_ref + 
                          kin_params.r32*dy_ref + 
                          kin_params.r33*dz_ref;

      v_tan_x = kin_params.vay * rz_target - kin_params.vaz * ry_target;
      v_tan_y = kin_params.vaz * rx_target - kin_params.vax * rz_target;
      v_tan_z = kin_params.vax * ry_target - kin_params.vay * rx_target;
    }

    wx[i] += (kin_params.vtx + v_tan_x) * p_val;
    wy[i] += (kin_params.vty + v_tan_y) * p_val;
    wz[i] += (kin_params.vtz + v_tan_z) * p_val;
  }
}

template <typename T>
__global__ void compute_cheap_dist_kernel(
                T * __restrict__ d,
                const T * __restrict__ x,
                const T * __restrict__ y,
                const T * __restrict__ z,
                const int lx, const int ly, const int lz,
                const int nel, const int local_iters,
                int* __restrict__ nchange)
{
    int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nel) return;

    int iter = 1;
    bool element_changed_ever = false;
    bool changed_local = true;

    const int lxy = lx * ly;
    const int lxyz = lx * ly * lz;
    const int e_offset = e * lxyz;

    while (changed_local && iter <= local_iters) {
        changed_local = false;

        // Loop over GLL nodes in this element
        for (int k = 0; k < lz; ++k) {
            for (int j = 0; j < ly; ++j) {
                for (int i = 0; i < lx; ++i) {
                    int idx1 = i + j * lx + k * lxy + e_offset;
                    T x1 = x[idx1];
                    T y1 = y[idx1];
                    T z1 = z[idx1];
                    T d1 = d[idx1];

                    int i0 = max(0, i - 1);
                    int i1 = min(lx - 1, i + 1);
                    int j0 = max(0, j - 1);
                    int j1 = min(ly - 1, j + 1);
                    int k0 = max(0, k - 1);
                    int k1 = min(lz - 1, k + 1);

                    // Neighbor check
                    for (int kk = k0; kk <= k1; ++kk) {
                       for (int jj = j0; jj <= j1; ++jj) {
                          for (int ii = i0; ii <= i1; ++ii) {
                             if (ii == i && jj == j && kk == k) continue;

                             int idx2 = ii + jj * lx + kk * lxy + e_offset;
                             T x2 = x[idx2];
                             T y2 = y[idx2];
                             T z2 = z[idx2];
                             T d2 = d[idx2];

                             T dist = sqrt((x1 - x2)*(x1 - x2) +
                                           (y1 - y2)*(y1 - y2) +
                                           (z1 - z2)*(z1 - z2));
                             T dtmp = d2 + dist;

                             if (dtmp < d1) {
                                 d1 = dtmp;
                                 d[idx1] = d1; // Update locally
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
        atomicAdd(nchange, 1);
    }
}

#endif // __COMMON_ALE_KINEMATICS_KERNEL_H__