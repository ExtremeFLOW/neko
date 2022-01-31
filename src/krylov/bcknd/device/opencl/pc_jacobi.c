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

#include "jacobi_kernel.cl.h"

/** 
 * Fortran wrapper for device pc jacobi
 */
void opencl_jacobi_update(void *d,
			  void *dxt, void *dyt, void *dzt,
			  void *G11, void *G22, void *G33,
			  void *G12, void *G13, void *G23,
			  int *nel, int *lx) {
  cl_int err;
   
  if (jacobi_program == NULL)
    opencl_kernel_jit(jacobi_kernel, (cl_program *) &jacobi_program);

  const int nb = (((*nel) * (*lx) * (*lx) * (*lx)) + 256 - 1) / 256; 
  const size_t global_item_size = 256 * nb;				       
  const size_t local_item_size = 256;				       
    
#define STR(X) #X
#define CASE(LX)                                                               \
  case LX:                                                                     \
    {                                                                          \
      cl_kernel kernel = clCreateKernel(jacobi_program,			       \
					STR(jacobi_kernel_lx##LX), &err);      \
      CL_CHECK(err);                                                           \
                                                                               \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &d));	       \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &dxt));      \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dyt));      \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dzt));      \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &G11));      \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &G22));      \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &G33));      \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &G12));      \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &G13));      \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &G23));      \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(int), nel));		       \
                                                                               \
                                                                               \
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,        \
				      kernel, 1, NULL, &global_item_size,      \
				      &local_item_size, 0, NULL, NULL));       \
    }                                                                          \
    break

  switch(*lx) {
    CASE(2);
    CASE(3);
    CASE(4);
    CASE(5);
    CASE(6);
    CASE(7);
    CASE(8);
    CASE(9);
    CASE(10);
    CASE(11);
    CASE(12);
    CASE(13);
    CASE(14);
    CASE(15);
    CASE(16);
  }
}
