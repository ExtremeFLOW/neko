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
#include <device/cuda/check.h>
#include <lpt/bcknd/device/cuda/lpt_wall_collision_kernel.h>

extern "C" {

  /** Fortran wrapper for the elastic particle-wall collision kernel */
  void cuda_lpt_wall_collide(
      void *xf, void *fr, void *fi, void *zg,
      void *el_ptr, void *el_facets, void *el_radius,
      int *lx, int *nelv, void *el_list, void *rst,
      void *x_old, void *y_old, void *z_old,
      void *x, void *y, void *z,
      void *d, void *u, void *v, void *w,
      void *u_lag, void *v_lag, void *w_lag,
      void *u_laglag, void *v_laglag, void *w_laglag,
      void *acc_xlag, void *acc_ylag, void *acc_zlag,
      void *acc_xlaglag, void *acc_ylaglag, void *acc_zlaglag,
      void *u_old, void *v_old, void *w_old,
      void *acc_x, void *acc_y, void *acc_z,
      int *lag_len, real *e_n, int *n, void *stats, cudaStream_t strm) {

    if (*n <= 0) {
      return;
    }

    const dim3 nthrds(128, 1, 1);
    const dim3 nblcks(((*n) + 128 - 1) / 128, 1, 1);

    lpt_wall_collide_kernel<real>
      <<<nblcks, nthrds, 0, strm>>>(
          (real *) xf, (real *) fr, (int *) fi, (real *) zg,
          (int *) el_ptr, (int *) el_facets, (real *) el_radius,
          *lx, *nelv, (int *) el_list, (real *) rst,
          (real *) x_old, (real *) y_old, (real *) z_old,
          (real *) x, (real *) y, (real *) z,
          (real *) d, (real *) u, (real *) v, (real *) w,
          (real *) u_lag, (real *) v_lag, (real *) w_lag,
          (real *) u_laglag, (real *) v_laglag, (real *) w_laglag,
          (real *) acc_xlag, (real *) acc_ylag, (real *) acc_zlag,
          (real *) acc_xlaglag, (real *) acc_ylaglag, (real *) acc_zlaglag,
          (real *) u_old, (real *) v_old, (real *) w_old,
          (real *) acc_x, (real *) acc_y, (real *) acc_z,
          *lag_len, *e_n, *n, (unsigned int *) stats);
    CUDA_CHECK(cudaGetLastError());
  }

}
