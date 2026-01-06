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
#include "compressible_ops_update_kernel.cl.h"

void opencl_update_uvw(void *u, void *v, void *w,
                       void *m_x, void *m_y, void *m_z,
                       void *rho, int n) {

  cl_int err;
  
  if (compressible_ops_update_program == NULL)
    opencl_kernel_jit(compressible_ops_update_kernel,
                      (cl_program *) &compressible_ops_update_program);

  cl_kernel kernel =
    clCreateKernel((cl_program) compressible_ops_update_program, 
                   "update_uvw_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m_x));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &m_y));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &m_z));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &rho));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), &n));

  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
} 

void opencl_update_mxyz_p_ruvw(void *m_x, void *m_y, void *m_z,
                               void *p, void *ruvw,
                               void *u, void *v, void *w, void *E,
                               void *rho, real gamma, int n) {
  cl_int err;
  
  if (compressible_ops_update_program == NULL)
    opencl_kernel_jit(compressible_ops_update_kernel,
                      (cl_program *) &compressible_ops_update_program);

  cl_kernel kernel =
    clCreateKernel((cl_program) compressible_ops_update_program, 
                   "update_mxyz_p_ruvw_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &m_x));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &m_y));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &m_z));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &p));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &ruvw));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &E));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &rho));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(real), &gamma));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(int), &n));

  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));  
}

void opencl_update_e(void *E, void *p, void *ruvw, real gamma, int n) {

  cl_int err;
  
  if (compressible_ops_update_program == NULL)
    opencl_kernel_jit(compressible_ops_update_kernel,
                      (cl_program *) &compressible_ops_update_program);

  cl_kernel kernel =
    clCreateKernel((cl_program) compressible_ops_update_program, 
                   "update_e_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &E));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &p));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &ruvw));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(real), &gamma));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), &n));

  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));  
}
