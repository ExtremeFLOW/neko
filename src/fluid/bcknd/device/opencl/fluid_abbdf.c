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

#include "abbdf_kernel.cl.h"

void fluid_sumab_opencl(void *u, void *v, void *w,
                        void *uu, void *vv, void *ww,
                        void *ulag1, void *ulag2, void *vlag1,
                        void *vlag2, void *wlag1, void *wlag2,
                        real *ab1, real *ab2, real *ab3, int *nab, int *n) {
  cl_int err;
  
  if (abbdf_program == NULL)
    opencl_kernel_jit(abbdf_kernel, (cl_program *) &abbdf_program);

  cl_kernel kernel = clCreateKernel(abbdf_program, "sumab_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &uu));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &vv));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &ww));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &ulag1));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &ulag2));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &vlag1));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &vlag2));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &wlag1));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &wlag2));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), ab1));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(real), ab2));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(real), ab3));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(int), nab));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

void fluid_makeabf_opencl(void *ta1, void *ta2, void *ta3,
                          void *abx1, void *aby1, void *abz1, 
                          void *abx2, void *aby2, void *abz2,
                          void *bfx, void *bfy, void *bfz,
                          real *rho, real *ab1, real *ab2, real *ab3, int *n) {
  cl_int err;
  
  if (abbdf_program == NULL)
    opencl_kernel_jit(abbdf_kernel, (cl_program *) &abbdf_program);
  
  cl_kernel kernel = clCreateKernel(abbdf_program, "makeabf_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &ta1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &ta2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &ta3));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &abx1));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &aby1));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &abz1));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &abx2));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &aby2));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &abz2));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bfx));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &bfy));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &bfz));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), rho));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(real), ab1));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(real), ab2));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(real), ab3));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
  
}

void fluid_makebdf_opencl(void *ta1, void *ta2, void *ta3,
                          void *tb1, void *tb2, void *tb3,
                          void *ulag1, void *ulag2, void *vlag1,
                          void *vlag2, void *wlag1, void *wlag2, 
                          void *bfx, void *bfy, void *bfz,
                          void *u, void *v, void *w, void *B, 
                          real *rho, real *dt, real *bd2,
                          real *bd3, real *bd4, int *nbd, int *n) {
  cl_int err;
  
  if (abbdf_program == NULL)
    opencl_kernel_jit(abbdf_kernel, (cl_program *) &abbdf_program);

  cl_kernel kernel = clCreateKernel(abbdf_program, "makebdf_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &ta1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &ta2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &ta3));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &tb1));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &tb2));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &tb3));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &ulag1));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &ulag2));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &vlag1));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &vlag2));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &wlag1));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &wlag2));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &bfx));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &bfy));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &bfz));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), (void *) &B));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(real), rho));
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(real), dt));
  CL_CHECK(clSetKernelArg(kernel, 21, sizeof(real), bd2));
  CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), bd3));
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(real), bd4));
  CL_CHECK(clSetKernelArg(kernel, 24, sizeof(int), nbd));
  CL_CHECK(clSetKernelArg(kernel, 25, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  
}
