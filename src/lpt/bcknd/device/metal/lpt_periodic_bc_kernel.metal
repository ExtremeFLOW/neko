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
 * Metal compute kernels for LPT periodic boundary wrapping.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Periodic direction data, packed to stay within Metal's kernel
 * argument limit. Each triple is (x, y, z) of direction @a j.
 * @note The layout must match lpt_periodic_params in
 * lpt/bcknd/device/lpt_periodic_params.h.
 */
struct lpt_periodic_params {
  float dir[9];
  float pmin[3];
  float pmax[3];
  float shift[9];
  float len[3];
};

/** Rotate the (x, y) pair of entry @a i by @a theta */
static void lpt_rotate_xy(device float *x, device float *y,
                          const int i, const float theta) {
  const float x_old = x[i];
  const float y_old = y[i];
  const float cos_theta = cos(theta);
  const float sin_theta = sin(theta);

  x[i] = cos_theta * x_old - sin_theta * y_old;
  y[i] = sin_theta * x_old + cos_theta * y_old;
}

kernel void lpt_periodic_bc_wrap_translational_kernel(
    device float *x [[ buffer(0) ]],
    device float *y [[ buffer(1) ]],
    device float *z [[ buffer(2) ]],
    constant int &n [[ buffer(3) ]],
    constant int &n_periodic_dirs [[ buffer(4) ]],
    constant lpt_periodic_params &prm [[ buffer(5) ]],
    uint gid [[ thread_position_in_grid ]]) {
  if (gid >= (uint) n) return;

  const int i = (int) gid;
  const float tol = 1.0e-8f;

  float point0 = x[i];
  float point1 = y[i];
  float point2 = z[i];

  for (int j = 0; j < n_periodic_dirs; j++) {
    const float dir0 = prm.dir[3 * j];
    const float dir1 = prm.dir[3 * j + 1];
    const float dir2 = prm.dir[3 * j + 2];
    const float pmin = prm.pmin[j];
    const float pmax = prm.pmax[j];
    const float shift0 = prm.shift[3 * j];
    const float shift1 = prm.shift[3 * j + 1];
    const float shift2 = prm.shift[3 * j + 2];
    const float len = prm.len[j];

    float coord = point0 * dir0 + point1 * dir1 + point2 * dir2;

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

/**
 * Rotational periodic wrapping.
 * @note @a present_mask selects which of the optional velocity and
 * acceleration pairs are bound; unset pairs are never dereferenced.
 */
kernel void lpt_periodic_bc_wrap_rotational_kernel(
    device float *x [[ buffer(0) ]],
    device float *y [[ buffer(1) ]],
    device float *z [[ buffer(2) ]],
    constant int &n [[ buffer(3) ]],
    constant float &theta_min [[ buffer(4) ]],
    constant float &theta_max [[ buffer(5) ]],
    constant float &theta_len [[ buffer(6) ]],
    device float *u [[ buffer(7) ]],
    device float *v [[ buffer(8) ]],
    device float *u_lag [[ buffer(9) ]],
    device float *v_lag [[ buffer(10) ]],
    device float *u_laglag [[ buffer(11) ]],
    device float *v_laglag [[ buffer(12) ]],
    device float *acc_xlag [[ buffer(13) ]],
    device float *acc_ylag [[ buffer(14) ]],
    device float *acc_xlaglag [[ buffer(15) ]],
    device float *acc_ylaglag [[ buffer(16) ]],
    constant int &present_mask [[ buffer(17) ]],
    uint gid [[ thread_position_in_grid ]]) {
  if (gid >= (uint) n) return;

  const int i = (int) gid;
  const float tol = 1.0e-8f;
  const float pi = acos(-1.0f);
  const float two_pi = 2.0f * pi;

  const float radius = sqrt(x[i] * x[i] + y[i] * y[i]);
  const float theta_old = fmod(atan2(y[i], x[i]) + two_pi, two_pi);
  float theta = theta_old;

  while (theta < theta_min - tol) {
    theta += theta_len;
  }

  while (theta > theta_max + tol) {
    theta -= theta_len;
  }

  const float dtheta = theta - theta_old;
  x[i] = radius * cos(theta);
  y[i] = radius * sin(theta);
  if (fabs(dtheta) <= tol) {
    return;
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
