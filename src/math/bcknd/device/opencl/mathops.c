/*
 Copyright (c) 2021-2022, The Neko Authors
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

#include "mathops_kernel.cl.h"

/** Fortran wrapper for opchsign \f$ a = -a \f$ */
void opencl_opchsign(void *a1, void *a2, void *a3, int *gdim, int *n) {    
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opchsign_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), gdim));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				  NULL, &global_item_size, &local_item_size,
				  0, NULL, NULL));
}

/** Fortran wrapper for opcolv \f$ a = a * c \f$ */
void opencl_opcolv(void *a1, void *a2, void *a3, void *c, int *gdim, int *n) {
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opcolv_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), gdim));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				  NULL, &global_item_size, &local_item_size,
				  0, NULL, NULL));   
}

/** Fortran wrapper for opcolv3c \f$ a(i) = b(i) * c(i) * d \f$ */
void opencl_opcolv3c(void *a1, void *a2, void *a3,
		     void *b1, void *b2, void *b3,
		     void *c, real *d, int *gdim, int *n) {
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opcolv3c_kernel", &err);
  CL_CHECK(err);
    
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &b1));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &b2));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &b3));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(real), d));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), gdim));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				  NULL, &global_item_size, &local_item_size,
				  0, NULL, NULL));   
}

/** Fortran wrapper for opadd2cm \f$ a(i) = a + b(i) * c \f$ */
void opencl_opadd2cm(void *a1, void *a2, void *a3, 
		     void *b1, void *b2, void *b3,
		     real *c, int *gdim, int *n) {
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opadd2cm_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &b1));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &b2));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &b3));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), gdim));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				  NULL, &global_item_size, &local_item_size,
				  0, NULL, NULL));  
}

/** Fortran wrapper for opadd2col \f$ a(i) = a + b(i) * c(i) \f$ */
void opencl_opadd2col(void *a1, void *a2, void *a3, 
		      void *b1, void *b2, void *b3,
		      void *c, int *gdim, int *n) {
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opadd2col_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &b1));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &b2));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &b3));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), gdim));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				  NULL, &global_item_size, &local_item_size,
				  0, NULL, NULL));
}
