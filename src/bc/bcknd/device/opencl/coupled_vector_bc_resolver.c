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
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>

#include "coupled_vector_bc_resolver_kernel.cl.h"

void opencl_coupled_vector_bc_resolver_apply(void *mixed_msk, void *x, void *y,
                                             void *z, void *constraint_n,
                                             void *constraint_t1,
                                             void *constraint_t2, void *n,
                                             void *t1, void *t2, int *m,
                                             cl_command_queue cmd_queue) {
  if (*m <= 0)
    return;

  cl_int err;

  if (coupled_vector_bc_resolver_program == NULL)
    opencl_kernel_jit(coupled_vector_bc_resolver_kernel,
                      (cl_program *) &coupled_vector_bc_resolver_program);

  cl_kernel kernel = clCreateKernel(coupled_vector_bc_resolver_program,
                                    "coupled_vector_bc_resolver_apply_kernel",
                                    &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &mixed_msk));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &y));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &z));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &constraint_n));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &constraint_t1));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &constraint_t2));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &n));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &t1));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &t2));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(int), m));

  const int nb = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}
