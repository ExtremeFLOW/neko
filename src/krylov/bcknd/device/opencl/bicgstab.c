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
 * OpenCL host-side dispatch for the fused BiCGStab kernels.
 *
 * The reductions return the process local sums; the Fortran side combines
 * them across ranks. A two value reduction keeps the partial sums of the
 * second value at @a buf_h[nb + group], so a single kernel launch and a
 * single read back cover both.
 */

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>
#include <device/opencl/buffer.h>

#include "bicgstab_kernel.cl.h"

/** Partial sums, up to two values per work group */
static opencl_buffer_t redbuf = OPENCL_BUFFER_INIT;

/**
 * Fortran wrapper for the BiCGStab search direction update
 * \f$ p = r + \beta (p - \omega v) \f$
 */
void opencl_bicgstab_update_p(void *p, void *r, void *v, real *beta,
                              real *omega, int *n,
                              cl_command_queue cmd_queue) {
  cl_int err;

  if (*n <= 0) {
    return;
  }

  if (bicgstab_program == NULL)
    opencl_kernel_jit(bicgstab_kernel, (cl_program *) &bicgstab_program);

  cl_kernel kernel = clCreateKernel(bicgstab_program,
                                    "bicgstab_update_p_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &p));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &r));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(real), beta));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(real), omega));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}

/**
 * Fortran wrapper for a weighted inner product and squared norm,
 * \f$ (a^T M b, b^T M b) \f$. The process local sums are returned in @a res.
 */
void opencl_bicgstab_product_and_norm(void *a, void *b, void *mult,
                                      real_xp *res, int *n,
                                      cl_command_queue cmd_queue) {
  cl_int err;
  cl_event kern_wait;
  int i;

  res[0] = 0.0;
  res[1] = 0.0;
  if (*n <= 0) {
    return;
  }

  if (bicgstab_program == NULL)
    opencl_kernel_jit(bicgstab_kernel, (cl_program *) &bicgstab_program);

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  opencl_buffer_reserve(&redbuf, 2 * nb * sizeof(real_xp));

  cl_kernel kernel = clCreateKernel(bicgstab_program,
                                    "bicgstab_product_and_norm_kernel",
                                    &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &redbuf.dev));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), n));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), &nb));

  CL_CHECK(clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, &kern_wait));

  CL_CHECK(clEnqueueReadBuffer(cmd_queue, redbuf.dev, CL_TRUE, 0,
                               2 * nb * sizeof(real_xp),
                               ((real_xp *) redbuf.host), 1,
                               &kern_wait, NULL));

  for (i = 0; i < nb; i++) {
    res[0] += ((real_xp *) redbuf.host)[i];
    res[1] += ((real_xp *) redbuf.host)[nb + i];
  }

  CL_CHECK(clReleaseEvent(kern_wait));
  CL_CHECK(clReleaseKernel(kernel));
}

/**
 * Fortran wrapper for BiCGStab part 1, \f$ s = r - \alpha v \f$,
 * returning the process local \f$ s^T M s \f$
 */
real_xp opencl_bicgstab_part1(void *s, void *r, void *v, void *mult,
                              real *alpha, int *n,
                              cl_command_queue cmd_queue) {
  cl_int err;
  cl_event kern_wait;
  int i;

  if (*n <= 0) {
    return 0.0;
  }

  if (bicgstab_program == NULL)
    opencl_kernel_jit(bicgstab_kernel, (cl_program *) &bicgstab_program);

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  opencl_buffer_reserve(&redbuf, nb * sizeof(real_xp));

  cl_kernel kernel = clCreateKernel(bicgstab_program,
                                    "bicgstab_part1_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &s));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &r));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &redbuf.dev));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(real), alpha));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(int), n));

  CL_CHECK(clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, &kern_wait));

  CL_CHECK(clEnqueueReadBuffer(cmd_queue, redbuf.dev, CL_TRUE, 0,
                               nb * sizeof(real_xp),
                               ((real_xp *) redbuf.host), 1,
                               &kern_wait, NULL));

  real_xp res = 0.0;
  for (i = 0; i < nb; i++) {
    res += ((real_xp *) redbuf.host)[i];
  }

  CL_CHECK(clReleaseEvent(kern_wait));
  CL_CHECK(clReleaseKernel(kernel));

  return res;
}

/**
 * Fortran wrapper for BiCGStab part 2,
 * \f$ x = x + \alpha \hat{p} + \omega \hat{s} \f$ and
 * \f$ r = s - \omega t \f$. The process local \f$ r^T M r \f$ and
 * \f$ f^T M r \f$ are returned in @a res.
 */
void opencl_bicgstab_part2(void *x, void *r, void *p_hat, void *s_hat,
                           void *s, void *t, void *f, void *mult,
                           real *alpha, real *omega, real_xp *res, int *n,
                           cl_command_queue cmd_queue) {
  cl_int err;
  cl_event kern_wait;
  int i;

  res[0] = 0.0;
  res[1] = 0.0;
  if (*n <= 0) {
    return;
  }

  if (bicgstab_program == NULL)
    opencl_kernel_jit(bicgstab_kernel, (cl_program *) &bicgstab_program);

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  opencl_buffer_reserve(&redbuf, 2 * nb * sizeof(real_xp));

  cl_kernel kernel = clCreateKernel(bicgstab_program,
                                    "bicgstab_part2_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &r));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &p_hat));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &s_hat));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &s));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &t));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &f));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &redbuf.dev));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(real), alpha));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(real), omega));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(int), n));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(int), &nb));

  CL_CHECK(clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, &kern_wait));

  CL_CHECK(clEnqueueReadBuffer(cmd_queue, redbuf.dev, CL_TRUE, 0,
                               2 * nb * sizeof(real_xp),
                               ((real_xp *) redbuf.host), 1,
                               &kern_wait, NULL));

  for (i = 0; i < nb; i++) {
    res[0] += ((real_xp *) redbuf.host)[i];
    res[1] += ((real_xp *) redbuf.host)[nb + i];
  }

  CL_CHECK(clReleaseEvent(kern_wait));
  CL_CHECK(clReleaseKernel(kernel));
}
