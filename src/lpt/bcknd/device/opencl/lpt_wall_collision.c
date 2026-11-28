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

/** Launch configuration shared by the two kernels */
static void lpt_wall_enqueue(cl_kernel kernel, int n, void *strm) {
  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) strm, kernel,
                                  1, NULL, &global_item_size,
                                  &local_item_size, 0, NULL, NULL));
}

static cl_kernel lpt_wall_kernel(const char *name) {
  cl_int err;
  if (lpt_wall_collision_program == NULL) {
    opencl_kernel_jit(lpt_wall_collision_kernel,
                      (cl_program *) &lpt_wall_collision_program);
  }
  cl_kernel kernel = clCreateKernel(lpt_wall_collision_program, name, &err);
  CL_CHECK(err);
  return kernel;
}

/**
 * Fortran wrapper for the OpenCL collision search: corrects the positions
 * and stores hit counts and contact normals per particle
 */
void opencl_lpt_wall_collide(void *xf, void *fr, void *fi, void *zg,
                             void *el_ptr, void *el_facets, void *el_radius,
                             int *lx, int *nelv, void *el_list, void *rst,
                             void *x_old, void *y_old, void *z_old,
                             void *x, void *y, void *z, void *d,
                             void *hits, void *normals, void *stats,
                             int *n, real *e_n, int *has_rst, void *strm) {
  if (*n < 1) return;

  cl_kernel kernel = lpt_wall_kernel("lpt_wall_collide_kernel");

  void *bufs[19] = { xf, fr, fi, zg, el_ptr, el_facets, el_radius,
                     el_list, rst, x_old, y_old, z_old, x, y, z, d,
                     hits, normals, stats };
  /* rst is only read when has_rst is set, but needs a valid object */
  if (bufs[8] == NULL) bufs[8] = x;

  for (int i = 0; i < 19; i++) {
    CL_CHECK(clSetKernelArg(kernel, i, sizeof(cl_mem), (void *) &bufs[i]));
  }
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(int), lx));
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(int), nelv));
  CL_CHECK(clSetKernelArg(kernel, 21, sizeof(int), n));
  CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), e_n));
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(int), has_rst));

  lpt_wall_enqueue(kernel, *n, strm);
  CL_CHECK(clReleaseKernel(kernel));
}

/**
 * Fortran wrapper for the OpenCL reflection of the velocity, the step's
 * velocity and acceleration, and their histories
 */
void opencl_lpt_wall_reflect(void *hits, void *normals,
                             void *u, void *v, void *w,
                             void *u_old, void *v_old, void *w_old,
                             void *acc_x, void *acc_y, void *acc_z,
                             void *u_lag, void *v_lag, void *w_lag,
                             void *acc_xlag, void *acc_ylag, void *acc_zlag,
                             void *u_laglag, void *v_laglag, void *w_laglag,
                             void *acc_xlaglag, void *acc_ylaglag,
                             void *acc_zlaglag, int *n, int *lag_len,
                             real *e_n, void *strm) {
  if (*n < 1) return;

  cl_kernel kernel = lpt_wall_kernel("lpt_wall_reflect_kernel");

  void *bufs[23] = { hits, normals, u, v, w, u_old, v_old, w_old,
                     acc_x, acc_y, acc_z, u_lag, v_lag, w_lag,
                     acc_xlag, acc_ylag, acc_zlag, u_laglag, v_laglag,
                     w_laglag, acc_xlaglag, acc_ylaglag, acc_zlaglag };

  /* The lag arrays are only read when lag_len allows it, but every slot
     still needs a valid memory object */
  for (int i = 0; i < 23; i++) {
    if (bufs[i] == NULL) bufs[i] = u;
  }

  for (int i = 0; i < 23; i++) {
    CL_CHECK(clSetKernelArg(kernel, i, sizeof(cl_mem), (void *) &bufs[i]));
  }
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(int), n));
  CL_CHECK(clSetKernelArg(kernel, 24, sizeof(int), lag_len));
  CL_CHECK(clSetKernelArg(kernel, 25, sizeof(real), e_n));

  lpt_wall_enqueue(kernel, *n, strm);
  CL_CHECK(clReleaseKernel(kernel));
}
