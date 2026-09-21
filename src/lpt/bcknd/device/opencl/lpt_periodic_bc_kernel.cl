#ifndef __LPT_PERIODIC_BC_KERNEL_CL__
#define __LPT_PERIODIC_BC_KERNEL_CL__
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
 * Periodic direction data, packed to stay within the device backends'
 * kernel argument limits. Each triple is (x, y, z) of direction @a j.
 * @note The layout must match lpt_periodic_params in the Metal backend
 * and in the host wrappers.
 */
typedef struct {
  real dir[9];
  real pmin[3];
  real pmax[3];
  real shift[9];
  real len[3];
} lpt_periodic_params;

/** Rotate the (x, y) pair of entry @a i by @a theta */
inline void lpt_rotate_xy(__global real * __restrict__ x,
                          __global real * __restrict__ y,
                          const int i, const real theta) {
  const real x_old = x[i];
  const real y_old = y[i];
  const real cos_theta = cos(theta);
  const real sin_theta = sin(theta);

  x[i] = cos_theta * x_old - sin_theta * y_old;
  y[i] = sin_theta * x_old + cos_theta * y_old;
}

/** Device kernel for translational periodic wrapping */
__kernel void lpt_periodic_bc_wrap_translational_kernel(
    __global real * __restrict__ x,
    __global real * __restrict__ y,
    __global real * __restrict__ z,
    const int n,
    const int n_periodic_dirs,
    const lpt_periodic_params prm) {

  const real tol = (real) 1.0e-8;
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    real point0 = x[i];
    real point1 = y[i];
    real point2 = z[i];

    for (int j = 0; j < n_periodic_dirs; j++) {
      const real dir0 = prm.dir[3 * j];
      const real dir1 = prm.dir[3 * j + 1];
      const real dir2 = prm.dir[3 * j + 2];
      const real pmin = prm.pmin[j];
      const real pmax = prm.pmax[j];
      const real shift0 = prm.shift[3 * j];
      const real shift1 = prm.shift[3 * j + 1];
      const real shift2 = prm.shift[3 * j + 2];
      const real len = prm.len[j];

      real coord = point0 * dir0 + point1 * dir1 + point2 * dir2;

      while (coord < pmin - tol) {
        point0 += shift0;
        point1 += shift1;
        point2 += shift2;
        coord += len;
      }

      while (coord > pmax + tol) {
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

/**
 * Device kernel for rotational periodic wrapping.
 * @note @a present_mask selects which of the optional velocity and
 * acceleration pairs are bound; unset pairs are never dereferenced.
 */
__kernel void lpt_periodic_bc_wrap_rotational_kernel(
    __global real * __restrict__ x,
    __global real * __restrict__ y,
    __global real * __restrict__ z,
    const int n,
    const real theta_min,
    const real theta_max,
    const real theta_len,
    __global real * __restrict__ u,
    __global real * __restrict__ v,
    __global real * __restrict__ u_lag,
    __global real * __restrict__ v_lag,
    __global real * __restrict__ u_laglag,
    __global real * __restrict__ v_laglag,
    __global real * __restrict__ acc_xlag,
    __global real * __restrict__ acc_ylag,
    __global real * __restrict__ acc_xlaglag,
    __global real * __restrict__ acc_ylaglag,
    const int present_mask) {

  const real tol = (real) 1.0e-8;
  const real pi = acos((real) -1.0);
  const real two_pi = (real) 2.0 * pi;
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    const real radius = sqrt(x[i] * x[i] + y[i] * y[i]);
    const real theta_old = fmod(atan2(y[i], x[i]) + two_pi, two_pi);
    real theta = theta_old;

    while (theta < theta_min - tol) {
      theta += theta_len;
    }

    while (theta > theta_max + tol) {
      theta -= theta_len;
    }

    const real dtheta = theta - theta_old;
    x[i] = radius * cos(theta);
    y[i] = radius * sin(theta);
    if (fabs(dtheta) <= tol) {
      continue;
    }

    if (present_mask & 1) {
      lpt_rotate_xy(u, v, i, dtheta);
    }
    if (present_mask & 2) {
      lpt_rotate_xy(u_lag, v_lag, i, dtheta);
    }
    if (present_mask & 4) {
      lpt_rotate_xy(u_laglag, v_laglag, i, dtheta);
    }
    if (present_mask & 8) {
      lpt_rotate_xy(acc_xlag, acc_ylag, i, dtheta);
    }
    if (present_mask & 16) {
      lpt_rotate_xy(acc_xlaglag, acc_ylaglag, i, dtheta);
    }
  }
}

#endif // __LPT_PERIODIC_BC_KERNEL_CL__
