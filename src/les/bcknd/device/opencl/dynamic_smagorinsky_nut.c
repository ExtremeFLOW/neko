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

#include "dynamic_smagorinsky_nut_kernel.cl.h"

/** JIT compile the dynamic Smagorinsky program on first use */
static void dynamic_smagorinsky_jit() {
  if (dynamic_smagorinsky_nut_program == NULL) {
    opencl_kernel_jit(dynamic_smagorinsky_nut_kernel,
                      (cl_program *) &dynamic_smagorinsky_nut_program);
  }
}

/** Enqueue a one-dimensional kernel over @a n points */
static void dynamic_smagorinsky_run(cl_kernel kernel, const int n) {
  if (n > 0) {
    const int nb = (n + 256 - 1) / 256;
    const size_t global_item_size = 256 * nb;
    const size_t local_item_size = 256;

    CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                    1, NULL, &global_item_size,
                                    &local_item_size, 0, NULL, NULL));
  }

  CL_CHECK(clReleaseKernel(kernel));
}

/** Fortran wrapper for the OpenCL strain rate magnitude kernel */
void opencl_s_abs_compute(void *s_abs, void *s11, void *s22, void *s33,
                          void *s12, void *s13, void *s23, int *n) {
  cl_int err;

  dynamic_smagorinsky_jit();

  cl_kernel kernel = clCreateKernel(dynamic_smagorinsky_nut_program,
                                    "s_abs_compute_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &s_abs));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &s11));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &s22));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &s33));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &s12));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &s13));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &s23));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), n));

  dynamic_smagorinsky_run(kernel, *n);
}

/** Fortran wrapper for the OpenCL Leonard stress kernel, part 1 */
void opencl_lij_compute_part1(void *l11, void *l22, void *l33,
                              void *l12, void *l13, void *l23,
                              void *u, void *v, void *w,
                              void *fu, void *fv, void *fw,
                              void *fuu, void *fvv, void *fww,
                              void *fuv, void *fuw, void *fvw, int *n) {
  cl_int err;

  dynamic_smagorinsky_jit();

  cl_kernel kernel = clCreateKernel(dynamic_smagorinsky_nut_program,
                                    "lij_compute_part1_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &l11));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &l22));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &l33));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &l12));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &l13));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &l23));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &fu));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &fv));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &fw));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &fuu));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &fvv));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &fww));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &fuv));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &fuw));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &fvw));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(int), n));

  dynamic_smagorinsky_run(kernel, *n);
}

/** Fortran wrapper for the OpenCL Leonard stress kernel, part 2 */
void opencl_lij_compute_part2(void *l11, void *l22, void *l33,
                              void *l12, void *l13, void *l23,
                              void *fuu, void *fvv, void *fww,
                              void *fuv, void *fuw, void *fvw, int *n) {
  cl_int err;

  dynamic_smagorinsky_jit();

  cl_kernel kernel = clCreateKernel(dynamic_smagorinsky_nut_program,
                                    "lij_compute_part2_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &l11));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &l22));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &l33));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &l12));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &l13));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &l23));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &fuu));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &fvv));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &fww));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &fuv));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &fuw));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &fvw));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(int), n));

  dynamic_smagorinsky_run(kernel, *n);
}

/** Fortran wrapper for the OpenCL model tensor kernel, part 1 */
void opencl_mij_compute_part1(void *m11, void *m22, void *m33,
                              void *m12, void *m13, void *m23,
                              void *s_abs, void *s11, void *s22, void *s33,
                              void *s12, void *s13, void *s23,
                              void *fs_abs, void *fs11, void *fs22, void *fs33,
                              void *fs12, void *fs13, void *fs23,
                              void *fsabss11, void *fsabss22, void *fsabss33,
                              void *fsabss12, void *fsabss13, void *fsabss23,
                              real *delta_ratio2, int *n) {
  cl_int err;

  dynamic_smagorinsky_jit();

  cl_kernel kernel = clCreateKernel(dynamic_smagorinsky_nut_program,
                                    "mij_compute_part1_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &m11));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &m22));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &m33));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m12));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &m13));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &m23));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &s_abs));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &s11));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &s22));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &s33));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &s12));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &s13));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &s23));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &fs_abs));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &fs11));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &fs22));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &fs33));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &fs12));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), (void *) &fs13));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), (void *) &fs23));
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), (void *) &fsabss11));
  CL_CHECK(clSetKernelArg(kernel, 21, sizeof(cl_mem), (void *) &fsabss22));
  CL_CHECK(clSetKernelArg(kernel, 22, sizeof(cl_mem), (void *) &fsabss33));
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(cl_mem), (void *) &fsabss12));
  CL_CHECK(clSetKernelArg(kernel, 24, sizeof(cl_mem), (void *) &fsabss13));
  CL_CHECK(clSetKernelArg(kernel, 25, sizeof(cl_mem), (void *) &fsabss23));
  CL_CHECK(clSetKernelArg(kernel, 26, sizeof(real), delta_ratio2));
  CL_CHECK(clSetKernelArg(kernel, 27, sizeof(int), n));

  dynamic_smagorinsky_run(kernel, *n);
}

/** Fortran wrapper for the OpenCL model tensor and eddy viscosity kernel */
void opencl_mij_nut_compute_part2(void *m11, void *m22, void *m33,
                                  void *m12, void *m13, void *m23,
                                  void *l11, void *l22, void *l33,
                                  void *l12, void *l13, void *l23,
                                  void *fsabss11, void *fsabss22,
                                  void *fsabss33, void *fsabss12,
                                  void *fsabss13, void *fsabss23,
                                  void *num, void *den, void *c_dyn,
                                  void *delta, void *s_abs, void *nut,
                                  real *alpha, int *n) {
  cl_int err;

  dynamic_smagorinsky_jit();

  cl_kernel kernel = clCreateKernel(dynamic_smagorinsky_nut_program,
                                    "mij_nut_compute_part2_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &m11));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &m22));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &m33));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m12));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &m13));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &m23));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &l11));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &l22));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &l33));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &l12));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &l13));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &l23));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &fsabss11));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &fsabss22));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &fsabss33));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &fsabss12));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &fsabss13));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &fsabss23));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), (void *) &num));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), (void *) &den));
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), (void *) &c_dyn));
  CL_CHECK(clSetKernelArg(kernel, 21, sizeof(cl_mem), (void *) &delta));
  CL_CHECK(clSetKernelArg(kernel, 22, sizeof(cl_mem), (void *) &s_abs));
  CL_CHECK(clSetKernelArg(kernel, 23, sizeof(cl_mem), (void *) &nut));
  CL_CHECK(clSetKernelArg(kernel, 24, sizeof(real), alpha));
  CL_CHECK(clSetKernelArg(kernel, 25, sizeof(int), n));

  dynamic_smagorinsky_run(kernel, *n);
}
