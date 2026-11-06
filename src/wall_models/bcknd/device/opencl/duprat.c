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

#include "duprat_kernel.cl.h"

/**
 * Fortran wrapper for the OpenCL kernel of the Duprat wall model.
 * @param u_d The x component of the sampled velocity.
 * @param v_d The y component of the sampled velocity.
 * @param w_d The z component of the sampled velocity.
 * @param n_x_d The x component of the wall normals.
 * @param n_y_d The y component of the wall normals.
 * @param n_z_d The z component of the wall normals.
 * @param nu_d The kinematic viscosity at wall points.
 * @param rho_w_d The density at wall points.
 * @param h_d The wall-model sampling distances.
 * @param dpds_d The wall-tangential pressure gradient (positive if adverse).
 * @param tau_x_d The x component of the wall shear stress.
 * @param tau_y_d The y component of the wall shear stress.
 * @param tau_z_d The z component of the wall shear stress.
 * @param alpha_d The ratio u_tau^2 / u_tau_p^2, for diagnostics.
 * @param n_nodes The number of wall points.
 * @param kappa The von Karman coefficient.
 * @param beta The exponent of the pressure-gradient term.
 * @param A The damping constant.
 * @param tstep The current time-step.
 */
void opencl_duprat_compute(void *u_d, void *v_d, void *w_d, void *n_x_d,
                           void *n_y_d, void *n_z_d, void *nu_d, void *rho_w_d,
                           void *h_d, void *dpds_d, void *tau_x_d,
                           void *tau_y_d, void *tau_z_d, void *alpha_d,
                           int *n_nodes, real *kappa, real *beta, real *A,
                           int *tstep) {
  cl_int err;

  if (*n_nodes < 1)
    return;

  if (duprat_program == NULL) {
    opencl_kernel_jit(duprat_kernel, (cl_program *)&duprat_program);
  }

  cl_kernel kernel =
      clCreateKernel(duprat_program, "duprat_compute_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&u_d));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&v_d));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&w_d));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&n_x_d));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&n_y_d));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&n_z_d));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&nu_d));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&rho_w_d));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&h_d));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&dpds_d));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&tau_x_d));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *)&tau_y_d));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *)&tau_z_d));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *)&alpha_d));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(int), n_nodes));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(real), kappa));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(real), beta));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(real), A));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(int), tstep));

  const int nb = ((*n_nodes) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue)glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size, 0,
                                  NULL, NULL));

  CL_CHECK(clReleaseKernel(kernel));
}

/**
 * Fortran wrapper for the OpenCL pressure-gradient update of the Duprat wall
 * model.
 * @param dpx_d The x component of the sampled pressure gradient.
 * @param dpy_d The y component of the sampled pressure gradient.
 * @param dpz_d The z component of the sampled pressure gradient.
 * @param u_d The x component of the sampled velocity.
 * @param v_d The y component of the sampled velocity.
 * @param w_d The z component of the sampled velocity.
 * @param n_x_d The x component of the wall normals.
 * @param n_y_d The y component of the wall normals.
 * @param n_z_d The z component of the wall normals.
 * @param nu_d The kinematic viscosity at wall points.
 * @param rho_w_d The density at wall points.
 * @param h_d The wall-model sampling distances.
 * @param eps The filter weight of the current gradient.
 * @param dpds_d The filtered pressure gradient, updated in place.
 * @param n_nodes The number of wall points.
 */
void opencl_duprat_update_dpds(void *dpx_d, void *dpy_d, void *dpz_d, void *u_d,
                               void *v_d, void *w_d, void *n_x_d, void *n_y_d,
                               void *n_z_d, void *nu_d, void *rho_w_d,
                               void *h_d, real *eps, void *dpds_d,
                               int *n_nodes) {
  cl_int err;

  if (*n_nodes < 1)
    return;

  if (duprat_program == NULL) {
    opencl_kernel_jit(duprat_kernel, (cl_program *)&duprat_program);
  }

  cl_kernel kernel =
      clCreateKernel(duprat_program, "duprat_update_dpds_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&dpx_d));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&dpy_d));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&dpz_d));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&u_d));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&v_d));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&w_d));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&n_x_d));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&n_y_d));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&n_z_d));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&nu_d));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&rho_w_d));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *)&h_d));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(real), eps));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *)&dpds_d));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(int), n_nodes));

  const int nb = ((*n_nodes) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue)glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size, 0,
                                  NULL, NULL));

  CL_CHECK(clReleaseKernel(kernel));
}
