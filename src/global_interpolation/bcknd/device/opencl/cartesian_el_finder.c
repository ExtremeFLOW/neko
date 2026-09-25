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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <device/device_config.h>
#include <device/opencl/check.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include "cartesian_el_finder_kernel.cl.h"

static void cartesian_el_finder_program_init(int xp_bytes) {
  if (cartesian_el_finder_program != NULL) return;
  if ((sizeof(real) != 4 && sizeof(real) != 8) ||
      (xp_bytes != 4 && xp_bytes != 8)) abort();

  const char *point_type = sizeof(real) == 4 ? "float" : "double";
  const char *xp_type = xp_bytes == 4 ? "float" : "double";
  const char *fp64 = (sizeof(real) == 8 || xp_bytes == 8) ?
    "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n" : "";
  char prefix[160];
  snprintf(prefix, sizeof(prefix), "%s#define POINT_T %s\n#define XP_T %s\n",
           fp64, point_type, xp_type);

  const size_t length = strlen(prefix) + strlen(cartesian_el_finder_kernel) + 1;
  char *source = (char *) malloc(length);
  if (source == NULL) abort();
  snprintf(source, length, "%s%s", prefix, cartesian_el_finder_kernel);
  opencl_kernel_jit(source, (cl_program *) &cartesian_el_finder_program);
  free(source);
}

void opencl_cartesian_el_finder_count(
    void *points, void *el_map_offset, void *point_box, void *n_el_cands,
    const void *min_x, const void *min_y, const void *min_z,
    const void *x_res, const void *y_res, const void *z_res,
    const int *n_boxes, const int *n_points, const int *xp_bytes) {
  if (*n_points <= 0) return;
  cartesian_el_finder_program_init(*xp_bytes);

  cl_int err;
  cl_kernel kernel = clCreateKernel(cartesian_el_finder_program,
                                    "cartesian_el_finder_count_kernel", &err);
  CL_CHECK(err);
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), &points));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), &el_map_offset));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), &point_box));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), &n_el_cands));
  CL_CHECK(clSetKernelArg(kernel, 4, *xp_bytes, min_x));
  CL_CHECK(clSetKernelArg(kernel, 5, *xp_bytes, min_y));
  CL_CHECK(clSetKernelArg(kernel, 6, *xp_bytes, min_z));
  CL_CHECK(clSetKernelArg(kernel, 7, *xp_bytes, x_res));
  CL_CHECK(clSetKernelArg(kernel, 8, *xp_bytes, y_res));
  CL_CHECK(clSetKernelArg(kernel, 9, *xp_bytes, z_res));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(int), n_boxes));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(int), n_points));

  const size_t local_size = 256;
  const size_t global_size = (((size_t)*n_points + local_size - 1) /
                              local_size) * local_size;
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                  1, NULL, &global_size, &local_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}

void opencl_cartesian_el_finder_fill(
    void *point_box, void *candidate_offsets, void *el_map_offset,
    void *el_map_data, void *candidate_array, const int *n_points) {
  if (*n_points <= 0) return;

  cl_int err;
  cl_kernel kernel = clCreateKernel(cartesian_el_finder_program,
                                    "cartesian_el_finder_fill_kernel", &err);
  CL_CHECK(err);
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), &point_box));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), &candidate_offsets));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), &el_map_offset));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), &el_map_data));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), &candidate_array));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n_points));

  const size_t local_size = 256;
  const size_t global_size = (((size_t)*n_points + local_size - 1) /
                              local_size) * local_size;
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                  1, NULL, &global_size, &local_size,
                                  0, NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}
