/*
 Copyright (c) 2024-2026, The Neko Authors
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
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>

#include "idw_kernel.cl.h"

/**
 * Fortran wrapper for the atomic-free IDW gather kernel.
 */
void opencl_idw_gather_one_sided(void *fu, void *fv, void *fw,
                                 void *fu_ib, void *fv_ib, void *fw_ib,
                                 void *fum_ib, void *fvm_ib, void *fwm_ib,
                                 void *x, void *y, void *z, void *ds,
                                 void *pmsk, void *w, void *wm,
                                 void *lpx, void *lpy, void *lpz,
                                 void *active_el, void *el_off, void *el_lag,
                                 int *n_active, int *lx3, real *dt,
                                 real *rmax, real *pwr, real *eps, real *wtol,
                                 cl_command_queue cmd_queue) {
  cl_int err;

  const int total = (*n_active) * (*lx3);
  if (total < 1)
    return;

  if (idw_program == NULL)
    opencl_kernel_jit(idw_kernel, (cl_program *) &idw_program);

  cl_kernel kernel = clCreateKernel(idw_program, "idw_gather_one_sided", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &fu));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &fv));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &fw));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &fu_ib));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &fv_ib));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &fw_ib));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &fum_ib));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &fvm_ib));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &fwm_ib));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &y));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &z));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &ds));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &pmsk));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &wm));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &lpx));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &lpy));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), (void *) &lpz));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), (void *) &active_el));
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), (void *) &el_off));
  CL_CHECK(clSetKernelArg(kernel, 21, sizeof(cl_mem), (void *) &el_lag));
  CL_CHECK(clSetKernelArg(kernel, 22, sizeof(int), n_active));
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(int), lx3));
  CL_CHECK(clSetKernelArg(kernel, 24, sizeof(real), dt));
  CL_CHECK(clSetKernelArg(kernel, 25, sizeof(real), rmax));
  CL_CHECK(clSetKernelArg(kernel, 26, sizeof(real), pwr));
  CL_CHECK(clSetKernelArg(kernel, 27, sizeof(real), eps));
  CL_CHECK(clSetKernelArg(kernel, 28, sizeof(real), wtol));

  const int nb = (total + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL,
                                  &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}
