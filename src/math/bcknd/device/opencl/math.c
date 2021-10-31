#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>

#include "math_kernel.cl.h"

/** Fortran wrapper for copy
 * Copy a vector \f$ a = b \f$
 */
void opencl_copy(void *a, void *b, int *n) {
  cl_int err = clEnqueueCopyBuffer((cl_command_queue) glb_cmd_queue,
				   a, b, 0, 0, (*n) * sizeof(real),
				   0, NULL, NULL);
}

/** Fortran wrapper for rzero
 * Zero a real vector
 */
void opencl_rzero(void *a, int *n) {
  real zero = 0;
  cl_int err = clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
  				   a, &zero, sizeof(real), 0,
				   (*n) * sizeof(real), 0, NULL, NULL);
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

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(real), c1);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
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

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(real), c1);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
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

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
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

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
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

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
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

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
  
  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
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

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
}

/**
 * Fortran wrapper glsc3
 * Weighted inner product \f$ a^T b c \f$
 */
real opencl_glsc3(void *a, void *b, void *c, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
    
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
    
  real * buf = (real *) malloc(nb * sizeof(real));

  cl_kernel kernel = clCreateKernel(math_program, "glsc3_kernel", &err);
  
  cl_mem buf_d = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE,
				nb * sizeof(real), NULL, &err);
  
  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &c);
  err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &buf_d);
  err = clSetKernelArg(kernel, 4, sizeof(int), n);
  
  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);

  err = clEnqueueReadBuffer((cl_command_queue) glb_cmd_queue, buf_d, CL_TRUE, 0,
			    nb * sizeof(real), buf, 0, NULL, NULL);
    
  real res = 0.0;
  for (int i = 0; i < nb; i++) {
    res += buf[i];
  }
  
  free(buf);
  err = clReleaseMemObject(buf_d);
  
  return res;
}

/**
 * Fortran wrapper glsc2
 * Weighted inner product \f$ a^T b c \f$
 */
real opencl_glsc2(void *a, void *b, int *n) {
  cl_int err;
  
  if (math_program == NULL)
    opencl_kernel_jit(math_kernel, (cl_program *) &math_program);
    
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
    
  real * buf = (real *) malloc(nb * sizeof(real));

  cl_kernel kernel = clCreateKernel(math_program, "glsc2_kernel", &err);
  
  cl_mem buf_d = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE,
				nb * sizeof(real), NULL, &err);
  
  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &buf_d);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);

  err = clEnqueueReadBuffer((cl_command_queue) glb_cmd_queue, buf_d, CL_TRUE, 0,
			    nb * sizeof(real), buf, 0, NULL, NULL);
    
  real res = 0.0;
  for (int i = 0; i < nb; i++) {
    res += buf[i];
  }
  
  free(buf);
  err = clReleaseMemObject(buf_d);
  
  return res;  
}
