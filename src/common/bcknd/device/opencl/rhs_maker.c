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

#include "rhs_maker_kernel.cl.h"

void rhs_maker_sumab_opencl(void *u, void *v, void *w,
                            void *uu, void *vv, void *ww,
                            void *ulag1, void *ulag2, void *vlag1,
                            void *vlag2, void *wlag1, void *wlag2,
                            real *ext1, real *ext2, real *ext3, int *nab, int *n) {
  cl_int err;
  
  if (rhs_maker_program == NULL)
    opencl_kernel_jit(rhs_maker_kernel, (cl_program *) &rhs_maker_program);

  cl_kernel kernel = clCreateKernel(rhs_maker_program, "sumab_kernel", &err);
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
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), ext1));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(real), ext2));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(real), ext3));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(int), nab));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

void rhs_maker_ext_opencl(void *abx1, void *aby1, void *abz1, 
                          void *abx2, void *aby2, void *abz2,
                          void *bfx, void *bfy, void *bfz,
                          real *rho, real *ext1, real *ext2, real *ext3, int *n) {
  cl_int err;
  
  if (rhs_maker_program == NULL)
    opencl_kernel_jit(rhs_maker_kernel, (cl_program *) &rhs_maker_program);
  
  cl_kernel kernel = clCreateKernel(rhs_maker_program, "makeext_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &abx1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &aby1));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &abz1));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &abx2));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &aby2));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &abz2));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &bfx));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &bfy));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &bfz));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(real), rho));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(real), ext1));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(real), ext2));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), ext3));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
  
}

void scalar_rhs_maker_ext_opencl(void *fs_lag, void *fs_laglag, void *fs,
                                 real *rho, real *ext1, real *ext2,
                                 real *ext3, int *n) {
  cl_int err;
  
  if (rhs_maker_program == NULL)
    opencl_kernel_jit(rhs_maker_kernel, (cl_program *) &rhs_maker_program);
  
  cl_kernel kernel = clCreateKernel(rhs_maker_program, "scalar_makeext_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &fs_lag));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &fs_laglag));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &fs));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(real), rho));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(real), ext1));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(real), ext2));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), ext3));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
  
}

void rhs_maker_bdf_opencl(void *ulag1, void *ulag2, void *vlag1,
                          void *vlag2, void *wlag1, void *wlag2, 
                          void *bfx, void *bfy, void *bfz,
                          void *u, void *v, void *w, void *B, 
                          real *rho, real *dt, real *bd2,
                          real *bd3, real *bd4, int *nbd, int *n) {
  cl_int err;
  
  if (rhs_maker_program == NULL)
    opencl_kernel_jit(rhs_maker_kernel, (cl_program *) &rhs_maker_program);

  cl_kernel kernel = clCreateKernel(rhs_maker_program, "makebdf_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &ulag1));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &ulag2));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &vlag1));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &vlag2));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &wlag1));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &wlag2));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &bfx));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &bfy));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &bfz));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &B));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(real), rho));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(real), dt));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(real), bd2));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(real), bd3));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(real), bd4));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(int), nbd));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  
}

void scalar_rhs_maker_bdf_opencl(void *s_lag, void *s_laglag, void *fs,
                                 void *s, void *B, real *rho, real *dt,
                                 real *bd2, real *bd3, real *bd4,
                                 int *nbd, int *n) {
  cl_int err;
  
  if (rhs_maker_program == NULL)
    opencl_kernel_jit(rhs_maker_kernel, (cl_program *) &rhs_maker_program);

  cl_kernel kernel = clCreateKernel(rhs_maker_program, "scalar_makebdf_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &s_lag));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &s_laglag));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &fs));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &s));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &B));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(real), rho));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(real), dt));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(real), bd2));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(real), bd3));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(real), bd4));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(int), nbd));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  
}
