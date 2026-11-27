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

#include "gradient_jump_penalty_kernel.cl.h"

#define GJP_NTHREADS 256

/** Fortran wrapper for the OpenCL facet value gather kernel */
void opencl_pick_facet_value_hex(void *b, void *a, int *nx, int *nel) {
  cl_int err;

  if (*nel < 1) return;

  if (gradient_jump_penalty_program == NULL) {
    opencl_kernel_jit(gradient_jump_penalty_kernel,
                      (cl_program *) &gradient_jump_penalty_program);
  }

  cl_kernel kernel = clCreateKernel(gradient_jump_penalty_program,
                                    "pick_facet_value_hex_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), nx));

  const size_t global_item_size = GJP_NTHREADS * (*nel);
  const size_t local_item_size = GJP_NTHREADS;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                  1, NULL, &global_item_size,
                                  &local_item_size, 0, NULL, NULL));

  CL_CHECK(clReleaseKernel(kernel));
}

/** Fortran wrapper for the OpenCL gradient jump penalty finalize kernel */
void opencl_gradient_jump_penalty_finalize(void *penalty_d,
                                           void *penalty_facet_d,
                                           void *dphidxi_d,
                                           int *nx, int *nel) {
  cl_int err;

  if (*nel < 1) return;

  if (gradient_jump_penalty_program == NULL) {
    opencl_kernel_jit(gradient_jump_penalty_kernel,
                      (cl_program *) &gradient_jump_penalty_program);
  }

  cl_kernel kernel =
    clCreateKernel(gradient_jump_penalty_program,
                   "gradient_jump_penalty_finalize_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &penalty_d));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem),
                          (void *) &penalty_facet_d));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dphidxi_d));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), nx));

  const size_t global_item_size = GJP_NTHREADS * (*nel);
  const size_t local_item_size = GJP_NTHREADS;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                  1, NULL, &global_item_size,
                                  &local_item_size, 0, NULL, NULL));

  CL_CHECK(clReleaseKernel(kernel));
}
