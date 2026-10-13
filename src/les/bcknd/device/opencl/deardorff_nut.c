/*
 Copyright (c) 2026, The Neko Authors
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

#include <device/device_config.h>
#include <device/opencl/check.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>

#include "deardorff_nut_kernel.cl.h"

/** Fortran wrapper for the OpenCL Deardorff eddy viscosity kernel */
void opencl_deardorff_nut_compute(void *TKE,
                                  void *dTdx, void *dTdy, void *dTdz,
                                  void *a11, void *a12, void *a13,
                                  void *a21, void *a22, void *a23,
                                  void *a31, void *a32, void *a33,
                                  void *delta, void *nut,
                                  void *temperature_alphat,
                                  void *TKE_alphat, void *TKE_source,
                                  real *c_k, real *T0,
                                  real *g1, real *g2, real *g3, real *eps,
                                  int *n) {
  cl_int err;

  if (deardorff_nut_program == NULL) {
    opencl_kernel_jit(deardorff_nut_kernel,
                      (cl_program *) &deardorff_nut_program);
  }

  cl_kernel kernel = clCreateKernel(deardorff_nut_program,
                                    "deardorff_nut_compute_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &TKE));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &dTdx));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dTdy));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dTdz));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &a11));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &a12));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &a13));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &a21));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &a22));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &a23));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &a31));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &a32));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &a33));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &delta));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &nut));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem),
                          (void *) &temperature_alphat));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &TKE_alphat));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &TKE_source));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(real), c_k));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(real), T0));
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(real), g1));
  CL_CHECK(clSetKernelArg(kernel, 21, sizeof(real), g2));
  CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), g3));
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(real), eps));
  CL_CHECK(clSetKernelArg(kernel, 24, sizeof(int), n));

  if (*n > 0) {
    const int nb = ((*n) + 256 - 1) / 256;
    const size_t global_item_size = 256 * nb;
    const size_t local_item_size = 256;

    CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                    1, NULL, &global_item_size,
                                    &local_item_size, 0, NULL, NULL));
  }

  CL_CHECK(clReleaseKernel(kernel));
}
