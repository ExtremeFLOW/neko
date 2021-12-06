#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>

#include "blasius_kernel.cl.h"

/** 
 * Fortran wrapper for device dirichlet apply vector
 */
void opencl_blasius_apply_vector(void *msk,
				 void *x, void *y, void *z,
				 void *bla_x, void *bla_y, void *bla_z,
				 int *m) {
  cl_int err;
  
  if (blasius_program == NULL)
    opencl_kernel_jit(blasius_kernel, (cl_program *) &blasius_program);
  
  cl_kernel kernel = clCreateKernel(blasius_program,
				    "blasius_apply_vector_kernel", &err);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &msk);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &x);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &y);
  err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &z);
  err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &bla_x);
  err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &bla_y);
  err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &bla_z);
  err = clSetKernelArg(kernel, 7, sizeof(int), m);
  
  const int nb = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);
}
