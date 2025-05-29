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

#include "euler_res_kernel.cl.h"

void euler_res_part_visc_opencl(void *rhs_u, void *Binv, void *lap_sol,
                                void *h, real *c_avisc, int *n) {
  cl_int err;
  
  if (euler_res_program == NULL)
    opencl_kernel_jit(euler_res_kernel, (cl_program *) &euler_res_program);

  cl_kernel kernel = clCreateKernel(euler_res_program,
                                    "euler_res_part_visc_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &rhs_u));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &Binv));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &lap_sol));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &h));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(real), c_avisc));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

void euler_res_part_mx_flux_opencl(void *f_x, void *f_y, void *f_z,
                                   void *m_x, void *m_y, void *m_z,
                                   void *rho_field, void *p, int *n) {
  cl_int err;
  
  if (euler_res_program == NULL)
    opencl_kernel_jit(euler_res_kernel, (cl_program *) &euler_res_program);
  
  cl_kernel kernel = clCreateKernel(euler_res_program,
                                    "euler_res_part_mx_flux_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &f_x));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &f_y));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &f_z));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m_x));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &m_y));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &m_z));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &rho_field));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &p));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

void euler_res_part_my_flux_opencl(void *f_x, void *f_y, void *f_z,
                                   void *m_x, void *m_y, void *m_z,
                                   void *rho_field, void *p, int *n) {
  cl_int err;
  
  if (euler_res_program == NULL)
    opencl_kernel_jit(euler_res_kernel, (cl_program *) &euler_res_program);
  
  cl_kernel kernel = clCreateKernel(euler_res_program,
                                    "euler_res_part_my_flux_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &f_x));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &f_y));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &f_z));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m_x));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &m_y));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &m_z));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &rho_field));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &p));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

void euler_res_part_mz_flux_opencl(void *f_x, void *f_y, void *f_z,
                                   void *m_x, void *m_y, void *m_z,
                                   void *rho_field, void *p, int *n) {
  cl_int err;
  
  if (euler_res_program == NULL)
    opencl_kernel_jit(euler_res_kernel, (cl_program *) &euler_res_program);
  
  cl_kernel kernel = clCreateKernel(euler_res_program,
                                    "euler_res_part_mz_flux_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &f_x));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &f_y));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &f_z));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m_x));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &m_y));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &m_z));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &rho_field));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &p));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

void euler_res_part_E_flux_opencl(void *f_x, void *f_y, void *f_z,
                                  void *m_x, void *m_y, void *m_z,
                                  void *rho_field, void *p, void * E, int *n) {
  cl_int err;
  
  if (euler_res_program == NULL)
    opencl_kernel_jit(euler_res_kernel, (cl_program *) &euler_res_program);
  
  cl_kernel kernel = clCreateKernel(euler_res_program,
                                    "euler_res_part_E_flux_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &f_x));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &f_y));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &f_z));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m_x));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &m_y));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &m_z));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &rho_field));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &E));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &p));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

void euler_res_part_coef_mult_opencl(void *rhs_rho, void *rhs_m_x,
                                     void *rhs_m_y, void *rhs_m_z,
                                     void *rhs_E, void *mult, int *n) {
  cl_int err;
  
  if (euler_res_program == NULL)
    opencl_kernel_jit(euler_res_kernel, (cl_program *) &euler_res_program);
  
  cl_kernel kernel = clCreateKernel(euler_res_program,
                                    "euler_res_part_coef_mult_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &rhs_rho));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &rhs_m_x));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &rhs_m_y));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &rhs_m_z));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &rhs_E));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

void euler_res_part_rk_sum_opencl(void *rho, void *m_x, void *m_y, void *m_z,
                                  void *E, void *k_rho_i, void *k_m_x_i,
                                  void *k_m_y_i, void *k_m_z_i, void *k_E_i,
                                  real *dt, real *c, int *n) {
  cl_int err;
  
  if (euler_res_program == NULL)
    opencl_kernel_jit(euler_res_kernel, (cl_program *) &euler_res_program);
  
  cl_kernel kernel = clCreateKernel(euler_res_program,
                                    "euler_res_part_rk_sum_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &rho));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &m_x));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &m_y));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m_z));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &E));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &k_rho_i));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &k_m_x_i));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &k_m_y_i));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &k_m_z_i));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &k_E_i));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(real), dt));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}
