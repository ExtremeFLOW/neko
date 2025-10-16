/*
 Copyright (c) 2025, The Neko Authors
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

#include <stdio.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>

#include "neumann_kernel.cl.h"

/** 
 * Fortran wrapper for device neumann apply scalar
 */
void opencl_neumann_apply_scalar(void *msk, void *facet,
                                 void *x, void *flux,
                                 void *area, int *lx, int *m,
                                 cl_command_queue cmd_queue) {

  cl_int err;
  
  if (neumann_program == NULL)
    opencl_kernel_jit(neumann_kernel, (cl_program *) &neumann_program);
  
  cl_kernel kernel = clCreateKernel(neumann_program,
                                  "neumann_apply_scalar_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &msk));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &facet));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &flux));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &area));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), lx));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(int), m));
  
  const int nb = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, NULL));
} 
