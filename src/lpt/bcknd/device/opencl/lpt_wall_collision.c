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

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <device/device_config.h>
#include <device/opencl/check.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>

#include "lpt_wall_collision_kernel.cl.h"

/** Fortran wrapper for the OpenCL elastic particle-wall collision kernel */
void opencl_lpt_handle_elastic_wall_collisions(
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
    int *lx, int *ly, int *lz, int *lag_len, void *strm) {
  cl_int err;

  if (*n < 1) return;

  if (lpt_wall_collision_program == NULL) {
    opencl_kernel_jit(lpt_wall_collision_kernel,
                      (cl_program *) &lpt_wall_collision_program);
  }

  cl_kernel kernel =
    clCreateKernel(lpt_wall_collision_program,
                   "lpt_handle_elastic_wall_collisions_kernel", &err);
  CL_CHECK(err);

  void *bufs[36] = { wall_facet_mask, el_list,
                     x_old, y_old, z_old,
                     x, y, z,
                     d, u, v, w,
                     u_lag, v_lag, w_lag,
                     u_laglag, v_laglag, w_laglag,
                     acc_xlag, acc_ylag, acc_zlag,
                     acc_xlaglag, acc_ylaglag, acc_zlaglag,
                     u_old, v_old, w_old,
                     acc_x, acc_y, acc_z,
                     dm_x, dm_y, dm_z,
                     nx, ny, nz };

  /* The lag arrays are only read when lag_len allows it, but every slot
     still needs a valid memory object */
  for (int i = 0; i < 36; i++) {
    if (bufs[i] == NULL) bufs[i] = x;
  }

  for (int i = 0; i < 36; i++) {
    CL_CHECK(clSetKernelArg(kernel, i, sizeof(cl_mem), (void *) &bufs[i]));
  }
  CL_CHECK(clSetKernelArg(kernel, 36, sizeof(int), n));
  CL_CHECK(clSetKernelArg(kernel, 37, sizeof(int), gdim));
  CL_CHECK(clSetKernelArg(kernel, 38, sizeof(int), nelv));
  CL_CHECK(clSetKernelArg(kernel, 39, sizeof(int), lx));
  CL_CHECK(clSetKernelArg(kernel, 40, sizeof(int), ly));
  CL_CHECK(clSetKernelArg(kernel, 41, sizeof(int), lz));
  CL_CHECK(clSetKernelArg(kernel, 42, sizeof(int), lag_len));

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) strm, kernel,
                                  1, NULL, &global_item_size,
                                  &local_item_size, 0, NULL, NULL));

  CL_CHECK(clReleaseKernel(kernel));
}
