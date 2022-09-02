/*
 Copyright (c) 2021-2022, The Neko Authors
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

#include "math_kernel.cl.h"

/** Fortran wrapper for copy
 * Copy a vector \f$ a = b \f$
 */
void opencl_copy(void *a, void *b, int *n) {
  CL_CHECK(clEnqueueCopyBuffer((cl_command_queue) glb_cmd_queue,
                               b, a, 0, 0, (*n) * sizeof(real),
                               0, NULL, NULL));
}

/** Fortran wrapper for rzero
 * Zero a real vector
 */
void opencl_rzero(void *a, int *n) {
  cl_event wait_kern;
  real zero = 0.0;

  CL_CHECK(clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
                               a, &zero, sizeof(real), 0,
                               (*n) * sizeof(real), 0, NULL, &wait_kern));
  CL_CHECK(clWaitForEvents(1, &wait_kern));
}

/** Fortran wrapper for rone
 * Set all elements to one
 */
void opencl_rone(void *a, int *n) {
  cl_event wait_kern;
  real one = 1.0;
  
  CL_CHECK(clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
                               a, &one, sizeof(real), 0,
                               (*n) * sizeof(real), 0, NULL, &wait_kern));
  CL_CHECK(clWaitForEvents(1, &wait_kern));
}

/** Fortran wrapper for cmult2
 * Multiplication by constant c \f$ a = c \cdot b \f$
 */
void opencl_cmult2(void *a, void *b, real *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
  
  cl_kernel kernel = clCreateKernel(math_program, "cmult2_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}


/** Fortran wrapper for cmult
 * Multiplication by constant c \f$ a = c \cdot a \f$
 */
void opencl_cmult(void *a, real *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
  
  cl_kernel kernel = clCreateKernel(math_program, "cmult_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}

/** Fortran wrapper for cadd
 * Add a scalar to vector \f$ a = \sum a_i + s \f$
 */
void opencl_cadd(void *a, real *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
  
  cl_kernel kernel = clCreateKernel(math_program, "cadd_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}

/** Fortran wrapper for cfill
 * Fill all elements to a constant c \f$ a = c  \f$
 */
void opencl_cfill(void *a, real *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
  
  cl_kernel kernel = clCreateKernel(math_program, "cfill_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}

/**
 * Fortran wrapper for add2
 * Vector addition \f$ a = a + b \f$
 */
void opencl_add2(void *a, void *b, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
  
  cl_kernel kernel = clCreateKernel(math_program, "add2_kernel", &err);
  CL_CHECK(err);
    
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}

/**
 * Fortran wrapper for add2s1
 * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
 * (multiplication on first argument) 
 */
void opencl_add2s1(void *a, void *b, real *c1, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
  
  cl_kernel kernel = clCreateKernel(math_program, "add2s1_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), c1));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}

/**
 * Fortran wrapper for add2s2
 * Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
 * (multiplication on second argument) 
 */
void opencl_add2s2(void *a, void *b, real *c1, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "add2s2_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), c1));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for add2s2
 * Vector addition with scalar multiplication 
 * \f$ x = x + c_1 p1 + c_2p2 + ... + c_jpj \f$
 * (multiplication on second argument) 
 */
void opencl_add2s2_many(void *x, void *p, void *alpha, int *j, int *n) {
  cl_int err;
  
  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "add2s2_many_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &p));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &alpha));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), j));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
  
}

/**
 * Fortran wrapper for addsqr2s2
 * Vector addition with scalar multiplication \f$ a = a + c_1 (b * b) \f$
 * (multiplication on second argument) 
 */
void opencl_addsqr2s2(void *a, void *b, real *c1, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "addsqr2s2_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), c1));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                               NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for add3s2
 * Vector addition with scalar multiplication \f$ a = c1 * b + c2 * c \f$
 */
void opencl_add3s2(void *a, void *b, void * c, real *c1, real *c2, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "add3s2_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(real), c1));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(real), c2));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for invcol1
 * Invert a vector \f$ a = 1 / a \f$
 */
void opencl_invcol1(void *a, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "invcol1_kernel", &err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(int), n));

  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}

/**
 * Fortran wrapper for invcol2
 * Vector division \f$ a = a / b \f$
 */
void opencl_invcol2(void *a, void *b, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "invcol2_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for col2
 * Vector multiplication with 2 vectors \f$ a = a \cdot b \f$
 */
void opencl_col2(void *a, void *b, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "col2_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for col3
 * Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
 */
void opencl_col3(void *a, void *b, void *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "col3_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}
  
/**
 * Fortran wrapper for subcol3
 * Vector multiplication with 3 vectors \f$ a = a - b \cdot c \f$
 */
void opencl_subcol3(void *a, void *b, void *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "subcol3_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for sub2
 * Vector subtraction \f$ a = a - b \f$
 */
void opencl_sub2(void *a, void *b, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "sub2_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for sub3
 * Vector subtraction \f$ a = b - c \f$
 */
void opencl_sub3(void *a, void *b, void *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "sub3_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for addcol3
 * \f$ a = a + b * c \f$
 */
void opencl_addcol3(void *a, void *b, void *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "addcol3_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/**
 * Fortran wrapper for addcol4
 * \f$ a = a + b * c * d \f$
 */
void opencl_addcol4(void *a, void *b, void *c, void *d, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);

  cl_kernel kernel = clCreateKernel(math_program, "addcol4_kernel", &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &d));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/** @todo cleanup this mess */
int red_s = 0;
real *bufred = NULL;
cl_mem bufred_d = NULL;

/**
 * Fortran wrapper glsc3
 * Weighted inner product \f$ a^T b c \f$
 */
real opencl_glsc3(void *a, void *b, void *c, int *n) {
  cl_int err;
  cl_event kern_wait;
  int i;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
    
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
    
  if ( nb > red_s){
    red_s = nb;
    if (bufred != NULL) {
      free(bufred);
      CL_CHECK(clReleaseMemObject(bufred_d));
    }
    bufred = (real *) malloc(nb * sizeof(real));
    bufred_d = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE,
                              nb * sizeof(real), NULL, &err);
    CL_CHECK(err);
  }
  
  cl_kernel kernel = clCreateKernel(math_program, "glsc3_kernel", &err);
  CL_CHECK(err);
      
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &bufred_d));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), n));
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, &kern_wait));
  
  CL_CHECK(clEnqueueReadBuffer((cl_command_queue) glb_cmd_queue, bufred_d,
                               CL_TRUE, 0, nb * sizeof(real), bufred, 1,
                               &kern_wait, NULL));
           
  real res = 0.0;
  for (i = 0; i < nb; i++) {
    res += bufred[i];
  }
  
  return res;
}

/**
 * Fortran wrapper for doing a reduction to an array
 * Weighted inner product \f$ w^T v(n,1:j) c \f$
 */
void opencl_glsc3_many(real *h, void * w, void *v, void *mult, int *j, int *n){ 
  int i, k;
  cl_int err;
  cl_event kern_wait;
  
  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
  
  int pow2 = 1;
  while(pow2 < (*j)){
    pow2 = 2*pow2;
  }

  const int nt = 256 / pow2;
  const int nb = ((*n) + nt - 1) / nt;
  const size_t local_item_size[2] = {nt, pow2};
  const size_t global_item_size[2] = {nb * nt, pow2};

  if((*j)*nb > red_s) {
    red_s = (*j)*nb;
    if (bufred != NULL) {
      free(bufred);    
      CL_CHECK(clReleaseMemObject(bufred_d));
    }
    bufred = (real *) malloc((*j)*nb * sizeof(real));
    bufred_d = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE,
                              (*j) * nb * sizeof(real), NULL, &err);
    CL_CHECK(err);
  }
  
  cl_kernel kernel = clCreateKernel(math_program, "glsc3_many_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &bufred_d));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), j));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 2,
                                  NULL, global_item_size, local_item_size,
                                  0, NULL, &kern_wait));
    
  CL_CHECK(clEnqueueReadBuffer((cl_command_queue) glb_cmd_queue,
                               bufred_d, CL_TRUE, 0, nb * sizeof(real),
                               bufred, 1, &kern_wait, NULL));

  for (k = 0; k < (*j); k++) {
    h[k] = 0.0;
  }
    
  for (i = 0; i < nb; i++) {
    for (k = 0; k < (*j); k++) {
        h[k] += bufred[i*(*j)+k];
    }
  }
}

/**
 * Fortran wrapper glsc2
 * Weighted inner product \f$ a^T b c \f$
 */
real opencl_glsc2(void *a, void *b, int *n) {
  cl_int err;
  cl_event kern_wait;
  int i;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
    
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
    
  real * buf = (real *) malloc(nb * sizeof(real));

  cl_kernel kernel = clCreateKernel(math_program, "glsc2_kernel", &err);
  CL_CHECK(err);
    
  cl_mem buf_d = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE,
                                nb * sizeof(real), NULL, &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &buf_d));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), n));
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, &kern_wait));

  
  CL_CHECK(clEnqueueReadBuffer((cl_command_queue) glb_cmd_queue, buf_d, CL_TRUE,
                               0, nb * sizeof(real), buf, 1, &kern_wait, NULL));
    
  real res = 0.0;
  for (i = 0; i < nb; i++) {
    res += buf[i];
  }

  free(buf);
  CL_CHECK(clReleaseMemObject(buf_d));
  
  return res;  
}

/**
 * Fortran wrapper glsum
 * Sum a vector of length n
 */
real opencl_glsum(void *a, int *n) {
  cl_int err;
  cl_event kern_wait;
  int i;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
    
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
    
  real * buf = (real *) malloc(nb * sizeof(real));

  cl_kernel kernel = clCreateKernel(math_program, "glsum_kernel", &err);
  CL_CHECK(err);
    
  cl_mem buf_d = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE,
                                nb * sizeof(real), NULL, &err);
  CL_CHECK(err);
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &buf_d));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), n));
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, &kern_wait));

  
  CL_CHECK(clEnqueueReadBuffer((cl_command_queue) glb_cmd_queue, buf_d, CL_TRUE,
                               0, nb * sizeof(real), buf, 1, &kern_wait, NULL));
    
  real res = 0.0;
  for (i = 0; i < nb; i++) {
    res += buf[i];
  }

  free(buf);
  CL_CHECK(clReleaseMemObject(buf_d));
  
  return res;  
}
