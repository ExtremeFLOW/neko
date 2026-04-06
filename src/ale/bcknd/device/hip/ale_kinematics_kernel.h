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
#endif // __COMMON_ALE_KINEMATICS_KERNEL_H__