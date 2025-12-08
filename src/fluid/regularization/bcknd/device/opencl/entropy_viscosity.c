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
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>

#include "entropy_viscosity_kernel.cl.h"

void opencl_entropy_visc_compute_residual(void *entropy_residual,
                                          void *S, void *S_lag1,
                                          void *S_lag2, void *S_lag3,
                                          real bdf1, real bdf2,
                                          real bdf3, real bdf4,
                                          real dt, int n) {
  cl_int err;

  if (entropy_viscosity_program == NULL)
    opencl_kernel_jit(entropy_viscosity_kernel,
                      (cl_program *) &entropy_viscosity_program);

  cl_kernel kernel = clCreateKernel(entropy_viscosity_program,
                                    "entropy_visc_compute_residual_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &entropy_residual));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &S));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &S_lag1));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &S_lag2));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &S_lag3));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(real), &bdf1));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(real), &bdf2));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(real), &bdf3));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(real), &bdf4));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(real), &dt));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(int), &n));

  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}

void opencl_entropy_visc_compute_viscosity(void *reg_coeff,
                                           void *entropy_residual,
                                           void *h,
                                           real c_entropy,
                                           real n_S,
                                           int n) {
  cl_int err;

  if (entropy_viscosity_program == NULL)
    opencl_kernel_jit(entropy_viscosity_kernel,
                      (cl_program *) &entropy_viscosity_program);

  cl_kernel kernel = clCreateKernel(entropy_viscosity_program,
                                    "entropy_visc_compute_viscosity_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &reg_coeff));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &entropy_residual));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &h));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(real), &c_entropy));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(real), &n_S));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), &n));

  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}

void opencl_entropy_visc_apply_element_max(void *reg_coeff,
                                           int lx,
                                           int nelv) {
  cl_int err;
  const int lx3 = lx * lx * lx;

  if (entropy_viscosity_program == NULL)
    opencl_kernel_jit(entropy_viscosity_kernel,
                      (cl_program *) &entropy_viscosity_program);

  cl_kernel kernel = clCreateKernel(entropy_viscosity_program,
                                    "entropy_visc_apply_element_max_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &reg_coeff));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(int), &lx3));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), &nelv));

  const int local_size = (lx3 < 256) ? lx3 : 256;
  const size_t global_item_size = local_size * nelv;
  const size_t local_item_size = local_size;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}

void opencl_entropy_visc_clamp_to_low_order(void *reg_coeff,
                                            void *h,
                                            void *max_wave_speed,
                                            real c_max,
                                            int n) {
  cl_int err;

  if (entropy_viscosity_program == NULL)
    opencl_kernel_jit(entropy_viscosity_kernel,
                      (cl_program *) &entropy_viscosity_program);

  cl_kernel kernel = clCreateKernel(entropy_viscosity_program,
                                    "entropy_visc_clamp_to_low_order_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &reg_coeff));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &h));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &max_wave_speed));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(real), &c_max));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), &n));

  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}

void opencl_entropy_visc_smooth_divide(void *reg_coeff,
                                       void *temp_field,
                                       void *mult_field,
                                       int n) {
  cl_int err;

  if (entropy_viscosity_program == NULL)
    opencl_kernel_jit(entropy_viscosity_kernel,
                      (cl_program *) &entropy_viscosity_program);

  cl_kernel kernel = clCreateKernel(entropy_viscosity_program,
                                    "entropy_visc_smooth_divide_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &reg_coeff));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &temp_field));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &mult_field));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), &n));

  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}

