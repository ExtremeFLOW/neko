/*
 Copyright (c) 2022, The Neko Authors
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

#include "schwarz_kernel.cl.h"

void opencl_schwarz_extrude(void *arr1, int * l1, real * f1,
                         void *arr2, int * l2, real * f2,
                         int * nx, int * nel) {  
  cl_int err;

  if (schwarz_program == NULL)
    opencl_kernel_jit(schwarz_kernel, (cl_program *) &schwarz_program);
  
  cl_kernel kernel = clCreateKernel(schwarz_program,
                                    "schwarz_extrude_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &arr1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(int), l1));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), f1));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &arr2));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), l1));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(real), f1));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(int), nel));
  
  const size_t global_item_size = 256 * (*nel);
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  
}

void opencl_schwarz_toext3d(void *a, void *b,int * nx, int * nel) {
  cl_int err;

  if (schwarz_program == NULL)
    opencl_kernel_jit(schwarz_kernel, (cl_program *) &schwarz_program);
  
  cl_kernel kernel = clCreateKernel(schwarz_program,
                                    "schwarz_toext3d_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), nx));
  
  const size_t global_item_size = 256 * (*nel);
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}

void opencl_schwarz_toreg3d(void *b, void *a,int * nx, int * nel) {
  cl_int err;

  if (schwarz_program == NULL)
    opencl_kernel_jit(schwarz_kernel, (cl_program *) &schwarz_program);
  
  cl_kernel kernel = clCreateKernel(schwarz_program,
                                    "schwarz_toreg3d_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), nx));
  
  const size_t global_item_size = 256 * (*nel);
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}


