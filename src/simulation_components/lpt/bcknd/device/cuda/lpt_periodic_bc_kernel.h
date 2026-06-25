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

#ifndef __LPT_PERIODIC_BC_KERNEL_H__
#define __LPT_PERIODIC_BC_KERNEL_H__

#include <math.h>

template<typename T>
__device__ inline void lpt_rotate_xy(T *x, T *y, const int i, const T theta) {
  const T x_old = x[i];
  const T y_old = y[i];
  const T cos_theta = cos(theta);
  const T sin_theta = sin(theta);

  x[i] = cos_theta * x_old - sin_theta * y_old;
  y[i] = sin_theta * x_old + cos_theta * y_old;
}

template<typename T>
__global__ void lpt_periodic_bc_wrap_translational_kernel(
    T * __restrict__ x,
    T * __restrict__ y,
    T * __restrict__ z,
    const int n,
    const int n_periodic_dirs,
    const T periodic_dir_x1,
    const T periodic_dir_y1,
    const T periodic_dir_z1,
    const T periodic_dir_x2,
    const T periodic_dir_y2,
    const T periodic_dir_z2,
    const T periodic_dir_x3,
    const T periodic_dir_y3,
    const T periodic_dir_z3,
    const T periodic_min1,
    const T periodic_min2,
    const T periodic_min3,
    const T periodic_max1,
    const T periodic_max2,
    const T periodic_max3,
    const T periodic_shift_x1,
    const T periodic_shift_y1,
    const T periodic_shift_z1,
    const T periodic_shift_x2,
    const T periodic_shift_y2,
    const T periodic_shift_z2,
    const T periodic_shift_x3,
    const T periodic_shift_y3,
    const T periodic_shift_z3,
    const T periodic_len1,
    const T periodic_len2,
    const T periodic_len3) {

  const T tol = 1.0e-8;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T point0 = x[i];
    T point1 = y[i];
    T point2 = z[i];

    for (int j = 0; j < n_periodic_dirs; j++) {
      T dir0 = periodic_dir_x1;
      T dir1 = periodic_dir_y1;
      T dir2 = periodic_dir_z1;
      T min = periodic_min1;
      T max = periodic_max1;
      T shift0 = periodic_shift_x1;
      T shift1 = periodic_shift_y1;
      T shift2 = periodic_shift_z1;
      T len = periodic_len1;

      switch (j) {
      case 1:
        dir0 = periodic_dir_x2;
        dir1 = periodic_dir_y2;
        dir2 = periodic_dir_z2;
        min = periodic_min2;
        max = periodic_max2;
        shift0 = periodic_shift_x2;
        shift1 = periodic_shift_y2;
        shift2 = periodic_shift_z2;
        len = periodic_len2;
        break;
      case 2:
        dir0 = periodic_dir_x3;
        dir1 = periodic_dir_y3;
        dir2 = periodic_dir_z3;
        min = periodic_min3;
        max = periodic_max3;
        shift0 = periodic_shift_x3;
        shift1 = periodic_shift_y3;
        shift2 = periodic_shift_z3;
        len = periodic_len3;
        break;
      default:
        break;
      }

      T coord = point0 * dir0 + point1 * dir1 + point2 * dir2;

      while (coord < min - tol) {
        point0 += shift0;
        point1 += shift1;
        point2 += shift2;
        coord += len;
      }

      while (coord > max + tol) {
        point0 -= shift0;
        point1 -= shift1;
        point2 -= shift2;
        coord -= len;
      }
    }

    x[i] = point0;
    y[i] = point1;
    z[i] = point2;
  }
}

template<typename T>
__global__ void lpt_periodic_bc_wrap_rotational_kernel(
    T * __restrict__ x,
    T * __restrict__ y,
    T * __restrict__ z,
    const int n,
    const T theta_min,
    const T theta_max,
    const T theta_len,
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
    T * __restrict__ acc_zlaglag) {

  const T tol = 1.0e-8;
  const T pi = acos((T) -1.0);
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const T radius = sqrt(x[i] * x[i] + y[i] * y[i]);
    const T theta_old = fmod(atan2(y[i], x[i]) + 2.0 * pi, 2.0 * pi);
    T theta = theta_old;

    while (theta < theta_min - tol) {
      theta += theta_len;
    }

    while (theta > theta_max + tol) {
      theta -= theta_len;
    }

    const T dtheta = theta - theta_old;
    x[i] = radius * cos(theta);
    y[i] = radius * sin(theta);
    if (fabs(dtheta) <= tol) {
      continue;
    }

    if (u && v) {
      lpt_rotate_xy(u, v, i, dtheta);
    }
    if (u_lag && v_lag) {
      lpt_rotate_xy(u_lag, v_lag, i, dtheta);
    }
    if (u_laglag && v_laglag) {
      lpt_rotate_xy(u_laglag, v_laglag, i, dtheta);
    }
    if (acc_xlag && acc_ylag) {
      lpt_rotate_xy(acc_xlag, acc_ylag, i, dtheta);
    }
    if (acc_xlaglag && acc_ylaglag) {
      lpt_rotate_xy(acc_xlaglag, acc_ylaglag, i, dtheta);
    }
  }
}

#endif
