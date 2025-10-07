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
#include <stdlib.h>
#include <math.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>

#include "prs_stress_res_kernel.cl.h"

void pnpn_prs_stress_res_part1_opencl(void *ta1, void *ta2, void *ta3,
                                      void *wa1, void *wa2, void *wa3,
                                      void *s11, void *s22, void *s33,
                                      void *s12, void *s13, void *s23,
                                      void *f_u, void *f_v, void *f_w,
                                      void *B, void *h1, void *rho, int *n) {

  cl_int err;

  if (pnpn_stress_res_program == NULL)
    opencl_kernel_jit(prs_stress_res_kernel,
                      (cl_program *) &pnpn_stress_res_program);

  cl_kernel kernel = clCreateKernel(pnpn_stress_res_program,
                                    "prs_stress_res_part1_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &ta1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &ta2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &ta3));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &wa1));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &wa2));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &wa3));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &s11));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &s22));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &s33));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &s12));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &s13));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &s23));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &f_u));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &f_v));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &f_w));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &B));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &rho));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(int), n));

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel(glb_cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, NULL));

}

void pnpn_prs_stress_res_part3_opencl(void *p_res, void *ta1, void *ta2,
                                      void *ta3, void *wa1, void *wa2,
                                      void *wa3, real *dtbd, int *n) {

  cl_int err;

  if (pnpn_stress_res_program == NULL)
    opencl_kernel_jit(prs_stress_res_kernel,
                      (cl_program *) &pnpn_stress_res_program);


  cl_kernel kernel = clCreateKernel(pnpn_stress_res_program,
                                    "prs_stress_res_part3_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &p_res));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &ta1));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &ta2));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ta3));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &wa1));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &wa2));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &wa3));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(real), dtbd));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), n));

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel(glb_cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, NULL));

}
