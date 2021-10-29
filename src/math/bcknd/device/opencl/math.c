#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "add2s1_kernel.cl.h"

/**
 * Fortran wrapper for add2s1
 * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
 * (multiplication on first argument) 
 */
void opencl_add2s1(cl_context ctx, cl_command_queue  cmd_queue,
		   cl_device_id device_id,
		   void *a, void *b, double *c1, int *n) {
  cl_int err;

  cl_program program =
    clCreateProgramWithSource(ctx, 1, &add2s1_kernel, NULL, &err);

  err = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
  
  cl_kernel kernel = clCreateKernel(program, "add2s1_kernel", &err);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &b);
  err = clSetKernelArg(kernel, 2, sizeof(double), c1);
  err = clSetKernelArg(kernel, 3, sizeof(int), n);
  
  const size_t global_item_size = CL_DEVICE_MAX_WORK_GROUP_SIZE;
  const size_t local_item_size = (((*n) + CL_DEVICE_MAX_WORK_GROUP_SIZE - 1) /
				  CL_DEVICE_MAX_WORK_GROUP_SIZE);

  err = clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL, 
			       &global_item_size, &local_item_size,
			       0, NULL, NULL);  
}
