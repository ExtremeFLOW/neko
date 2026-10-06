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

#include <math.h>

#include <device/device_config.h>
#include <device/opencl/check.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>

#include "vreman_nut_kernel.cl.h"

/** Fortran wrapper for the OpenCL Vreman eddy viscosity kernel */
void opencl_vreman_nut_compute(void *a11, void *a12, void *a13,
                               void *a21, void *a22, void *a23,
                               void *a31, void *a32, void *a33,
                               void *delta, void *nut, void *mult,
                               real *c, real *eps, int *n) {
  cl_int err;

  if (vreman_nut_program == NULL) {
    opencl_kernel_jit(vreman_nut_kernel, (cl_program *) &vreman_nut_program);
  }

  cl_kernel kernel = clCreateKernel(vreman_nut_program,
                                    "vreman_nut_compute_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a11));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a12));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a13));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &a21));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &a22));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &a23));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &a31));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &a32));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &a33));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &delta));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &nut));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(real), eps));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(int), n));

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

/**
 * Fortran wrapper for the OpenCL buoyancy-corrected Vreman kernel
 * @note @a g is the (unnormalised) gravity vector; the shear direction
 * uses the normalised vector, the buoyancy term the original one.
 */
void opencl_vreman_nut_compute_buoy(void *a11, void *a12, void *a13,
                                    void *a21, void *a22, void *a23,
                                    void *a31, void *a32, void *a33,
                                    void *delta, void *nut, void *mult,
                                    real *c, real *eps, int *n,
                                    void *dTdx, void *dTdy, void *dTdz,
                                    real *g, real *ri_c, real *ref_temp) {
  cl_int err;

  if (vreman_nut_program == NULL) {
    opencl_kernel_jit(vreman_nut_kernel, (cl_program *) &vreman_nut_program);
  }

  cl_kernel kernel = clCreateKernel(vreman_nut_program,
                                    "vreman_nut_compute_buoy_kernel", &err);
  CL_CHECK(err);

  const real g1 = g[0];
  const real g2 = g[1];
  const real g3 = g[2];
  const real gmag = sqrt(g1 * g1 + g2 * g2 + g3 * g3);
  const real n1 = g1 / gmag;
  const real n2 = g2 / gmag;
  const real n3 = g3 / gmag;

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a11));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a12));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a13));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &a21));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &a22));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &a23));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &a31));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &a32));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &a33));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &delta));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &nut));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(real), eps));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(int), n));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &dTdx));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &dTdy));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &dTdz));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(real), (void *) &n1));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(real), (void *) &n2));
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(real), (void *) &n3));
  CL_CHECK(clSetKernelArg(kernel, 21, sizeof(real), (void *) &g1));
  CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), (void *) &g2));
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(real), (void *) &g3));
  CL_CHECK(clSetKernelArg(kernel, 24, sizeof(real), ri_c));
  CL_CHECK(clSetKernelArg(kernel, 25, sizeof(real), ref_temp));

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
