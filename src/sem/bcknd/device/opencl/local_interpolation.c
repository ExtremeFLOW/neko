/*
 Copyright (c) 2022-2025, The Neko Authors
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>
#include <common/neko_log.h>

#include "local_interpolation_kernel.cl.h"

/**
 * Fortran wrapper for local interpolation
 */
void opencl_find_rst_legendre(void *rst,
                              void *pt_x, void* pt_y, void* pt_z,
                              void *x_hat, void *y_hat, void *z_hat,
                              void *resx, void *resy, void *resz,
                              int *lx, void *el_ids, int *n_pt, real *tol,
                              void *conv_pts) {
  cl_int err;
  if (find_rst_legendre_program == NULL)
    opencl_kernel_jit(local_interpolation_kernel,
                      (cl_program *) &find_rst_legendre_program);

  size_t global_kstep[2];
  size_t local_kstep[2];
  local_kstep[0] = 1;
  local_kstep[1] = 128;
  global_kstep[0] = (*n_pt);
  global_kstep[1] = 128;
  
#define STR(X) #X
#define CASE(LX)                                                               \
  case LX:                                                                     \
    {                                                                          \
      cl_kernel kernel =                                                       \
        clCreateKernel(find_rst_legendre_program,                              \
                       STR(find_rst_legendre_kernel_lx##LX), &err);            \
      CL_CHECK(err);                                                           \
                                                                               \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &rst));      \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &pt_x));     \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &pt_y));     \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &pt_z));     \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &x_hat));    \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &y_hat));    \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &z_hat));    \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &resx));     \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &resy));     \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &resz));     \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &el_ids));  \
      CL_CHECK(clSetKernelArg(kernel, 11, sizeof(int), n_pt));                 \
      CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), tol));                 \
      CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &conv_pts));\
                                                                               \
     CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,         \
                                     kernel, 2, NULL, global_kstep,            \
                                     local_kstep, 0, NULL, NULL));             \
                                                                               \
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
  default:
    {
      fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
  }
}
  
