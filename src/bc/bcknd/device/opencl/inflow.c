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

#include "inflow_kernel.cl.h"

/** 
 * Fortran wrapper for device inflow apply vector
 */
void opencl_inflow_apply_vector(void *msk, void *x, void *y,
				void *z, void *g, int *m) {

  cl_int err;
   
  if (inflow_program == NULL)
    opencl_kernel_jit(inflow_kernel, (cl_program *) &inflow_program);
  
  cl_kernel kernel = clCreateKernel(inflow_program,
				    "inflow_apply_vector_kernel", &err);
  CL_CHECK(err);

  const real gx = ((real *)g)[0];
  const real gy = ((real *)g)[1];
  const real gz = ((real *)g)[2];
  
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &msk));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &y));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &z));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(real), &gx));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(real), &gy));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(real), &gz));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), m));
  
  const int nb = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				  NULL, &global_item_size, &local_item_size,
				  0, NULL, NULL));

}
