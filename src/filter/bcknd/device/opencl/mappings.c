/*
 Copyright (c) 2021-2024, The Neko Authors
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
#include <stdio.h>
#include <stdlib.h>

#include "mapping_kernels.cl.h"

/** Fortran wrapper for 2nd order smooth_step
 *
 * Compute the smooth step function for a given array.
 * t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
 * return t^3 * (t * (6.0 * t - 15.0) + 10.0);
 */
void opencl_smooth_step(void* x, real* edge0, real* edge1, int* n) {
    cl_int err;

    if (mapping_program == NULL)
        opencl_kernel_jit(mapping_kernels, (cl_program*)mapping_program);

    cl_kernel kernel =
        clCreateKernel(mapping_program, "smooth_step_kernel", &err);
    CL_CHECK(err);

    CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&x));
    CL_CHECK(clSetKernelArg(kernel, 1, sizeof(real), edge0));
    CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), edge1));
    CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));

    const int    nb               = ((*n) + 256 - 1) / 256;
    const size_t global_item_size = 256 * nb;
    const size_t local_item_size  = 256;

    CL_CHECK(clEnqueueNDRangeKernel(
        (cl_command_queue)glb_cmd_queue, kernel, 1, NULL, &global_item_size,
        &local_item_size, 0, NULL, NULL));
}

/** Fortran wrapper for step function
 *
 * Compute the step function for a given array.
 * if (x < edge) return left;
 * else return right;
 *
 * @param x array to apply the step function
 * @param edge threshold value
 * @param left value to return if x < edge
 * @param right value to return if x >= edge
 * @param n size of the array
 */
void opencl_step_function(
    void* x, real* edge, real* left, real* right, int* n) {
    cl_int err;

    if (mapping_program == NULL)
        opencl_kernel_jit(mapping_kernels, (cl_program*)mapping_program);

    cl_kernel kernel =
        clCreateKernel(mapping_program, "step_function_kernel", &err);
    CL_CHECK(err);

    CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&x));
    CL_CHECK(clSetKernelArg(kernel, 1, sizeof(real), edge));
    CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), left));
    CL_CHECK(clSetKernelArg(kernel, 3, sizeof(real), right));
    CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), n));

    const int    nb               = ((*n) + 256 - 1) / 256;
    const size_t global_item_size = 256 * nb;
    const size_t local_item_size  = 256;

    CL_CHECK(clEnqueueNDRangeKernel(
        (cl_command_queue)glb_cmd_queue, kernel, 1, NULL, &global_item_size,
        &local_item_size, 0, NULL, NULL));
}

/** Fortran wrapper for the permeability mapping
 *
 * @param x array to apply the permeability mapping on
 * @param k_0 lower bound of the permeability
 * @param k_1 upper bound of the permeability
 * @param q parameter
 * @param n size of the array
 */
void opencl_permeability(void* x, real* k_0, real* k_1, real* q, int* n) {
    cl_int err;

    if (mapping_program == NULL)
        opencl_kernel_jit(mapping_kernels, (cl_program*)mapping_program);

    cl_kernel kernel =
        clCreateKernel(mapping_program, "permeability_kernel", &err);
    CL_CHECK(err);

    CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&x));
    CL_CHECK(clSetKernelArg(kernel, 1, sizeof(real), k_0));
    CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), k_1));
    CL_CHECK(clSetKernelArg(kernel, 3, sizeof(real), q));
    CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), n));

    const int    nb               = ((*n) + 256 - 1) / 256;
    const size_t global_item_size = 256 * nb;
    const size_t local_item_size  = 256;

    CL_CHECK(clEnqueueNDRangeKernel(
        (cl_command_queue)glb_cmd_queue, kernel, 1, NULL, &global_item_size,
        &local_item_size, 0, NULL, NULL));
}