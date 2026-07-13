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

#include <device/device_config.h>
#include <device/hip/check.h>
#include <lpt/bcknd/device/hip/lpt_wall_collision_kernel.h>

extern "C" {

  void hip_lpt_handle_elastic_wall_collisions(
      void *wall_facet_mask, void *el_list,
      void *x_old, void *y_old, void *z_old,
      void *x, void *y, void *z,
      void *d, void *u, void *v, void *w,
      void *u_lag, void *v_lag, void *w_lag,
      void *u_laglag, void *v_laglag, void *w_laglag,
      void *acc_xlag, void *acc_ylag, void *acc_zlag,
      void *acc_xlaglag, void *acc_ylaglag, void *acc_zlaglag,
      void *u_old, void *v_old, void *w_old,
      void *acc_x, void *acc_y, void *acc_z,
      void *dm_x, void *dm_y, void *dm_z,
      void *nx, void *ny, void *nz,
      int *n, int *gdim, int *nelv,
      int *lx, int *ly, int *lz, int *lag_len, hipStream_t strm) {

    if (*n <= 0) {
      return;
    }

    const dim3 nthrds(256, 1, 1);
    const dim3 nblcks(((*n) + 256 - 1) / 256, 1, 1);

    hipLaunchKernelGGL(
        HIP_KERNEL_NAME(lpt_handle_elastic_wall_collisions_kernel<real>),
        nblcks, nthrds, 0, strm,
        (int *) wall_facet_mask, (int *) el_list,
        (real *) x_old, (real *) y_old, (real *) z_old,
        (real *) x, (real *) y, (real *) z,
        (real *) d, (real *) u, (real *) v, (real *) w,
        (real *) u_lag, (real *) v_lag, (real *) w_lag,
        (real *) u_laglag, (real *) v_laglag, (real *) w_laglag,
        (real *) acc_xlag, (real *) acc_ylag, (real *) acc_zlag,
        (real *) acc_xlaglag, (real *) acc_ylaglag, (real *) acc_zlaglag,
        (real *) u_old, (real *) v_old, (real *) w_old,
        (real *) acc_x, (real *) acc_y, (real *) acc_z,
        (real *) dm_x, (real *) dm_y, (real *) dm_z,
        (real *) nx, (real *) ny, (real *) nz,
        *n, *gdim, *nelv, *lx, *ly, *lz, *lag_len);
    HIP_CHECK(hipGetLastError());
  }

}
