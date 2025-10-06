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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>
#include <common/neko_log.h>

#include "ax_helm_full_kernel.cl.h"

/**
 * Fortran wrapper for device OpenCL Ax helm full
 */
void opencl_ax_helm_stress_vector(void *au, void *av, void *aw,
                                  void *u, void *v, void *w,
                                  void *dx, void *dy, void *dz,
                                  void *dxt, void *dyt, void *dzt,
                                  void *h1,
                                  void *drdx, void *drdy, void *drdz,
                                  void *dsdx, void *dsdy, void *dsdz,
                                  void *dtdx, void *dtdy, void *dtdz,
                                  void *jacinv, void *w3, int *nelv, int *lx) {
  cl_int err;

  if (ax_helm_full_program == NULL)
    opencl_kernel_jit(ax_helm_full_kernel,
                      (cl_program *) &ax_helm_full_program);

  size_t global_kstep[2];
  size_t local_kstep[2];
  local_kstep[0] = (*lx);
  local_kstep[1] = (*lx);
  global_kstep[0] = (*nelv) * (*lx);
  global_kstep[1] = (*lx);

#define STR(X) #X
#define CASE(LX)                                                               \
  case LX:                                                                     \
    {                                                                          \
      cl_kernel kernel =                                                       \
        clCreateKernel(ax_helm_full_program,                                   \
                       STR(ax_helm_stress_kernel_full_lx##LX), &err);          \
      CL_CHECK(err);                                                           \
                                                                               \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &au));       \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &av));       \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &aw));       \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &u));        \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &v));        \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &w));        \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dx));       \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dy));       \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &dz));       \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &h1));       \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &drdx));    \
      CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &drdy));    \
      CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &drdz));    \
      CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &dsdx));    \
      CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &dsdy));    \
      CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &dsdz));    \
      CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &dtdx));    \
      CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &dtdy));    \
      CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), (void *) &dtdz));    \
      CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), (void *) &jacinv));  \
      CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), (void *) &w3));      \
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

  /**
   * Fortran wrapper for device OpenCL Ax helm full version part2
   */
void opencl_ax_helm_stress_vector_part2(void *au, void *av, void *aw,
                                        void *u, void *v, void *w,
                                        void *h2, void *B, int *n) {

  cl_int err;

  if (ax_helm_full_program == NULL)
    opencl_kernel_jit(ax_helm_full_kernel,
                      (cl_program *) &ax_helm_full_program);

  cl_kernel kernel = clCreateKernel(ax_helm_full_program,
                                    "ax_helm_stress_kernel_vector_part2", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &au));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &av));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &aw));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &h2));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &B));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), n));

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel(glb_cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, NULL));
 
}
