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

#include "most_kernel.cl.h"

/**
 * Fortran wrapper for the OpenCL most wall model kernel
 * @note @a tstep is unused; it is kept for interface compatibility
 * with the other wall models.
 */
void opencl_most_compute(void *u_d, void *v_d, void *w_d, void *temp_d,
                         void *temp_w_d,
                         void *n_x_d, void *n_y_d, void *n_z_d, void *h_d,
                         void *tau_x_d, void *tau_y_d, void *tau_z_d,
                         int *n_nodes, real *kappa, void *mu_w_d,
                         void *rho_w_d, real *g, real *Pr, real *z0,
                         real *z0h_in, int *bc_type, real *bc_value,
                         int *tstep,
                         void *Ri_b_diagn, void *L_ob_diagn, void *utau_diagn,
                         void *magu_diagn, void *ti_diagn, void *ts_diagn,
                         void *q_diagn) {
  cl_int err;

  if (*n_nodes < 1) return;
  if (*bc_type != 0 && *bc_type != 1) return;

  if (most_program == NULL) {
    opencl_kernel_jit(most_kernel, (cl_program *) &most_program);
  }

  cl_kernel kernel =
    clCreateKernel(most_program,
                   (*bc_type == 0) ? "most_compute_neumann_kernel"
                                   : "most_compute_dirichlet_kernel", &err);
  CL_CHECK(err);

  const real g1 = g[0];
  const real g2 = g[1];
  const real g3 = g[2];

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &u_d));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &v_d));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &w_d));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &temp_d));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &temp_w_d));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &h_d));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &n_x_d));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &n_y_d));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &n_z_d));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &tau_x_d));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &tau_y_d));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &tau_z_d));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(int), n_nodes));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(real), kappa));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &mu_w_d));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &rho_w_d));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(real), (void *) &g1));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(real), (void *) &g2));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(real), (void *) &g3));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(real), Pr));
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(real), z0));
  CL_CHECK(clSetKernelArg(kernel, 21, sizeof(real), z0h_in));
  CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), bc_value));
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(cl_mem), (void *) &Ri_b_diagn));
  CL_CHECK(clSetKernelArg(kernel, 24, sizeof(cl_mem), (void *) &L_ob_diagn));
  CL_CHECK(clSetKernelArg(kernel, 25, sizeof(cl_mem), (void *) &utau_diagn));
  CL_CHECK(clSetKernelArg(kernel, 26, sizeof(cl_mem), (void *) &magu_diagn));
  CL_CHECK(clSetKernelArg(kernel, 27, sizeof(cl_mem), (void *) &ti_diagn));
  CL_CHECK(clSetKernelArg(kernel, 28, sizeof(cl_mem), (void *) &ts_diagn));
  CL_CHECK(clSetKernelArg(kernel, 29, sizeof(cl_mem), (void *) &q_diagn));

  const int nb = ((*n_nodes) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                  1, NULL, &global_item_size,
                                  &local_item_size, 0, NULL, NULL));

  CL_CHECK(clReleaseKernel(kernel));
}
