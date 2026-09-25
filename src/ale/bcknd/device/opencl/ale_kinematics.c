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

#include "ale_kinematics_kernel.cl.h"

/**
 * Rigid body kinematics, passed by value.
 * @note The layout must match kinematics_params_t in the device kernels
 * and the bind(c) type in ale_routines_device.F90.
 */
typedef struct {
  real cx, cy, cz;
  real vtx, vty, vtz;
  real vax, vay, vaz;
  real px, py, pz;
  real r11, r12, r13;
  real r21, r22, r23;
  real r31, r32, r33;
} kinematics_params_t;

/** Fortran wrapper for the OpenCL mesh velocity kinematics kernel */
void add_kinematics_to_mesh_velocity_opencl(void *wx, void *wy, void *wz,
                                            void *x_ref, void *y_ref,
                                            void *z_ref, void *phi,
                                            void *x, void *y, void *z,
                                            kinematics_params_t kin_params,
                                            int n) {
  cl_int err;

  if (n < 1) return;

  if (ale_kinematics_program == NULL) {
    opencl_kernel_jit(ale_kinematics_kernel,
                      (cl_program *) &ale_kinematics_program);
  }

  cl_kernel kernel = clCreateKernel(ale_kinematics_program,
                                    "ale_add_kinematics_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(int), (void *) &n));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &wx));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &wy));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &wz));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &x_ref));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &y_ref));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &z_ref));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &phi));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &y));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &z));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(kinematics_params_t),
                          (void *) &kin_params));

  const int nb = (n + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                  1, NULL, &global_item_size,
                                  &local_item_size, 0, NULL, NULL));

  CL_CHECK(clReleaseKernel(kernel));
}

/** Fortran wrapper for the OpenCL cheap distance kernel */
void compute_cheap_dist_opencl(void *d_d, void *x_d, void *y_d, void *z_d,
                               int lx, int ly, int lz, int nel,
                               int local_iters, void *nchange_d) {
  cl_int err;

  if (nel < 1) return;

  if (ale_kinematics_program == NULL) {
    opencl_kernel_jit(ale_kinematics_kernel,
                      (cl_program *) &ale_kinematics_program);
  }

  cl_kernel kernel = clCreateKernel(ale_kinematics_program,
                                    "compute_cheap_dist_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &d_d));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &x_d));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &y_d));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &z_d));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), (void *) &lx));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), (void *) &ly));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(int), (void *) &lz));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), (void *) &nel));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), (void *) &local_iters));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &nchange_d));

  const int nb = (nel + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                  1, NULL, &global_item_size,
                                  &local_item_size, 0, NULL, NULL));

  CL_CHECK(clReleaseKernel(kernel));
}
