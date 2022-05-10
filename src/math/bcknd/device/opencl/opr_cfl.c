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

#include "cfl_kernel.cl.h"

/** 
 * Fortran wrapper for device OpenCL convective terms
 */
real opencl_cfl(real *dt, void *u, void *v, void *w,
		void *drdx, void *dsdx, void *dtdx,
		void *drdy, void *dsdy, void *dtdy,
		void *drdz, void *dsdz, void *dtdz,
		void *dr_inv, void *ds_inv, void *dt_inv,
		void *jacinv, int *nel, int *lx) {
  cl_int err;
  cl_event kern_wait;
  int i;
  if (cfl_program == NULL)
    opencl_kernel_jit(cfl_kernel, (cl_program *) &cfl_program);
  
  const size_t global_item_size = 256 * (*nel);
  const size_t local_item_size = 256;

  real * cfl = (real *) malloc((*nel) * sizeof(real));
  cl_mem cfl_d = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE,
                                (*nel) * sizeof(real), NULL, &err);
  CL_CHECK(err);

#define STR(X) #X
#define CASE(LX)                                                                \
  case LX:                                                                      \
    {                                                                           \
      cl_kernel kernel = clCreateKernel(cfl_program,                            \
                                        STR(cfl_kernel_lx##LX), &err);          \
      CL_CHECK(err);                                                            \
                                                                                \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(real), dt));                    \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u));         \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &v));         \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &w));         \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &drdx));      \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dsdx));      \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dtdx));      \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &drdy));      \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &dsdy));      \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &dtdy));      \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &drdz));     \
      CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &dsdz));     \
      CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &dtdz));     \
      CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &dr_inv));   \
      CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &ds_inv));   \
      CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &dt_inv));   \
      CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &jacinv));   \
      CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &cfl_d));    \
                                                                                \
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,         \
                                      kernel, 1, NULL, &global_item_size,       \
                                      &local_item_size, 0, NULL, &kern_wait));  \
    }                                                                           \
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
  }

  CL_CHECK(clEnqueueReadBuffer((cl_command_queue) glb_cmd_queue,
			       cfl_d, CL_TRUE, 0, (*nel) * sizeof(real),
			       cfl, 1, &kern_wait, NULL));
    
  real cfl_max = 0.0;
  for (i = 0; i < (*nel); i++) {
    cfl_max = fmax(cfl_max, cfl[i]);
  }

  free(cfl);
  CL_CHECK(clReleaseMemObject(cfl_d));

  return cfl_max;

} 
