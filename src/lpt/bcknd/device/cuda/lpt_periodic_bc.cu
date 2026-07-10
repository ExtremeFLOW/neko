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
#include <lpt/bcknd/device/cuda/lpt_periodic_bc_kernel.h>

extern "C" {

  void cuda_lpt_periodic_bc_wrap_translational(
      void *x, void *y, void *z, int *n, int *n_periodic_dirs,
      real *periodic_dir_x1, real *periodic_dir_y1, real *periodic_dir_z1,
      real *periodic_dir_x2, real *periodic_dir_y2, real *periodic_dir_z2,
      real *periodic_dir_x3, real *periodic_dir_y3, real *periodic_dir_z3,
      real *periodic_min1, real *periodic_min2, real *periodic_min3,
      real *periodic_max1, real *periodic_max2, real *periodic_max3,
      real *periodic_shift_x1, real *periodic_shift_y1,
      real *periodic_shift_z1, real *periodic_shift_x2,
      real *periodic_shift_y2, real *periodic_shift_z2,
      real *periodic_shift_x3, real *periodic_shift_y3,
      real *periodic_shift_z3, real *periodic_len1, real *periodic_len2,
      real *periodic_len3, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    lpt_periodic_bc_wrap_translational_kernel<real>
      <<<nblcks, nthrds, 0, strm>>>((real *) x, (real *) y, (real *) z,
                                    *n, *n_periodic_dirs, *periodic_dir_x1,
                                    *periodic_dir_y1, *periodic_dir_z1,
                                    *periodic_dir_x2, *periodic_dir_y2,
                                    *periodic_dir_z2, *periodic_dir_x3,
                                    *periodic_dir_y3, *periodic_dir_z3,
                                    *periodic_min1, *periodic_min2,
                                    *periodic_min3, *periodic_max1,
                                    *periodic_max2, *periodic_max3,
                                    *periodic_shift_x1, *periodic_shift_y1,
                                    *periodic_shift_z1, *periodic_shift_x2,
                                    *periodic_shift_y2, *periodic_shift_z2,
                                    *periodic_shift_x3, *periodic_shift_y3,
                                    *periodic_shift_z3, *periodic_len1,
                                    *periodic_len2, *periodic_len3);
    CUDA_CHECK(cudaGetLastError());
  }

  void cuda_lpt_periodic_bc_wrap_rotational(
      void *x, void *y, void *z, int *n,
      real *theta_min, real *theta_max, real *theta_len,
      void *u, void *v, void *w,
      void *u_lag, void *v_lag, void *w_lag,
      void *u_laglag, void *v_laglag, void *w_laglag,
      void *acc_xlag, void *acc_ylag, void *acc_zlag,
      void *acc_xlaglag, void *acc_ylaglag, void *acc_zlaglag,
      cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    lpt_periodic_bc_wrap_rotational_kernel<real>
      <<<nblcks, nthrds, 0, strm>>>((real *) x, (real *) y, (real *) z,
                                    *n, *theta_min, *theta_max, *theta_len,
                                    (real *) u, (real *) v, (real *) w,
                                    (real *) u_lag, (real *) v_lag,
                                    (real *) w_lag, (real *) u_laglag,
                                    (real *) v_laglag, (real *) w_laglag,
                                    (real *) acc_xlag, (real *) acc_ylag,
                                    (real *) acc_zlag,
                                    (real *) acc_xlaglag,
                                    (real *) acc_ylaglag,
                                    (real *) acc_zlaglag);
    CUDA_CHECK(cudaGetLastError());
  }

}
