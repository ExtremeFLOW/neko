#ifndef __COMMON_ALE_KINEMATICS_KERNEL_H__
#define __COMMON_ALE_KINEMATICS_KERNEL_H__

#include <cmath>
#include <algorithm>

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
                const T cx, const T cy, const T cz,
                const T vtx, const T vty, const T vtz,
                const T vax, const T vay, const T vaz,
                const T px, const T py, const T pz,
                const T r11, const T r12, const T r13,
                const T r21, const T r22, const T r23,
                const T r31, const T r32, const T r33) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const T p_val = phi[i];
    T v_tan_x, v_tan_y, v_tan_z;

    if (abs(p_val - 1.0) < 1e-6) {
      const T rx = x[i] - cx;
      const T ry = y[i] - cy;
      const T rz = z[i] - cz;

      v_tan_x = vay * rz - vaz * ry;
      v_tan_y = vaz * rx - vax * rz;
      v_tan_z = vax * ry - vay * rx;
    } 
    else {
      const T dx_ref = x_ref[i] - px;
      const T dy_ref = y_ref[i] - py;
      const T dz_ref = z_ref[i] - pz;

      const T rx_target = r11*dx_ref + r12*dy_ref + r13*dz_ref;
      const T ry_target = r21*dx_ref + r22*dy_ref + r23*dz_ref;
      const T rz_target = r31*dx_ref + r32*dy_ref + r33*dz_ref;

      v_tan_x = vay * rz_target - vaz * ry_target;
      v_tan_y = vaz * rx_target - vax * rz_target;
      v_tan_z = vax * ry_target - vay * rx_target;
    }

    wx[i] += (vtx + v_tan_x) * p_val;
    wy[i] += (vty + v_tan_y) * p_val;
    wz[i] += (vtz + v_tan_z) * p_val;
  }
}

template <typename T>
__global__ void compute_cheap_dist_kernel(
                T* __restrict__ d,
                const T* __restrict__ x,
                const T* __restrict__ y,
                const T* __restrict__ z,
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