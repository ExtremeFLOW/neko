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

#ifndef __LPT_WALL_COLLISION_KERNEL_H__
#define __LPT_WALL_COLLISION_KERNEL_H__

/*
 * Elastic collisions of particles with the high-order wall facets: one
 * thread per particle, the algorithm is shared with the other device
 * backends in lpt_wall_intersect.h (a port of lpt_wall_geometry.f90).
 */

#include <lpt/bcknd/device/lpt_wall_intersect.h>

template<typename T>
__global__ void lpt_wall_collide_kernel(
    const T * __restrict__ xf,
    const T * __restrict__ fr,
    const int * __restrict__ fi,
    const T * __restrict__ zg,
    const int * __restrict__ el_ptr,
    const int * __restrict__ el_facets,
    const T * __restrict__ el_radius,
    const int lx,
    const int nelv,
    const int * __restrict__ el_list,
    const T * __restrict__ rst,
    const T * __restrict__ x_old,
    const T * __restrict__ y_old,
    const T * __restrict__ z_old,
    T * __restrict__ x,
    T * __restrict__ y,
    T * __restrict__ z,
    const T * __restrict__ d,
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
    T * __restrict__ acc_zlaglag,
    T * __restrict__ u_old,
    T * __restrict__ v_old,
    T * __restrict__ w_old,
    T * __restrict__ acc_x,
    T * __restrict__ acc_y,
    T * __restrict__ acc_z,
    const int lag_len,
    const T e_n,
    const int n,
    unsigned int * __restrict__ stats) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  lpt_wall_geom<T> g;
  g.xf = xf;
  g.fr = fr;
  g.fi = fi;
  g.zg = zg;
  g.el_ptr = el_ptr;
  g.el_facets = el_facets;
  g.el_radius = el_radius;
  g.lx = lx;
  g.nelv = nelv;

  for (int p = idx; p < n; p += str) {
    lpt_wall_cnt cnt;
    for (int i = 0; i < LPT_NSTATS; i++) cnt.c[i] = 0;

    lpt_wall_collide_particle<T>(g, p, el_list, rst, x_old, y_old, z_old,
                                 x, y, z, d, u, v, w, u_lag, v_lag, w_lag,
                                 u_laglag, v_laglag, w_laglag, acc_xlag,
                                 acc_ylag, acc_zlag, acc_xlaglag,
                                 acc_ylaglag, acc_zlaglag, u_old, v_old,
                                 w_old, acc_x, acc_y, acc_z, lag_len, e_n,
                                 cnt);

    /* Most particles never come near a wall: only touch the global
       counters when something was counted beyond the particle itself. */
    for (int i = 0; i < LPT_NSTATS; i++) {
      if (cnt.c[i] == 0) continue;
      atomicAdd(&stats[i], cnt.c[i]);
    }
  }
}

#endif // __LPT_WALL_COLLISION_KERNEL_H__
