#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <device/device_config.h>

#include "math_kernel.cl.h"

cl_program math_program = NULL;

/**
 * JIT compile OpenCL math kernels
 */

void opencl_math_kernel_jit()
{
  cl_int err;

  math_program = clCreateProgramWithSource((cl_context) glb_ctx, 1,
					   &math_kernel, NULL, &err);
  if (err != CL_SUCCESS) printf("%d\n", err);
  err = clBuildProgram(math_program, 1, (cl_device_id *) &glb_device_id,
		       NULL, NULL, NULL);

  if (err == CL_BUILD_PROGRAM_FAILURE) {
    // Determine the size of the log
    size_t log_size;
    clGetProgramBuildInfo(math_program, glb_device_id,
			  CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    // Allocate memory for the log
    char *log = (char *) malloc(log_size);

    // Get the log
    clGetProgramBuildInfo(math_program, glb_device_id,
			  CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

    // Print the log
    printf("%s\n", log);

  }
  
}

/**
 * Fortran wrapper for add2s1
 * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
 * (multiplication on first argument) 
 */
void opencl_add2s1(void *a, void *b, real *c1, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_math_kernel_jit();
  
  cl_kernel kernel = clCreateKernel(math_program, "add2s1_kernel", &err);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &b);
  err = clSetKernelArg(kernel, 2, sizeof(real), c1);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  const size_t global_item_size = CL_DEVICE_MAX_WORK_GROUP_SIZE;
  const size_t local_item_size = (((*n) + CL_DEVICE_MAX_WORK_GROUP_SIZE - 1) /
				  CL_DEVICE_MAX_WORK_GROUP_SIZE);

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
    opencl_math_kernel_jit();

  cl_kernel kernel = clCreateKernel(math_program, "add2s2_kernel", &err);
  if (err != CL_SUCCESS) printf("%d\n", err);
  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &b);
  err = clSetKernelArg(kernel, 2, sizeof(real), c1);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  const size_t global_item_size = CL_DEVICE_MAX_WORK_GROUP_SIZE;
  const size_t local_item_size = (((*n) + CL_DEVICE_MAX_WORK_GROUP_SIZE - 1) /
				  CL_DEVICE_MAX_WORK_GROUP_SIZE);

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
  if (err != CL_SUCCESS) printf("%d\n", err);
}

/**
 * Fortran wrapper for col2
 * Vector multiplication with 2 vectors \f$ a = a \cdot b \f$
 */
void opencl_col2(void *a, void *b, int *n) {
  cl_int err;

  if (math_program == NULL)
    opencl_math_kernel_jit();

  cl_kernel kernel = clCreateKernel(math_program, "col2_kernel", &err);
  if (err != CL_SUCCESS) printf("%d\n", err);
  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &b);
  err = clSetKernelArg(kernel, 2, sizeof(int), n);
  
  const size_t global_item_size = CL_DEVICE_MAX_WORK_GROUP_SIZE;
  const size_t local_item_size = (((*n) + CL_DEVICE_MAX_WORK_GROUP_SIZE - 1) /
				  CL_DEVICE_MAX_WORK_GROUP_SIZE);

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
  if (err != CL_SUCCESS) printf("%d\n", err);
}
