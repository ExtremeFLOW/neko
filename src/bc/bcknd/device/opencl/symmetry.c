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

#include "symmetry_kernel.cl.h"

#define MAX(a, b) (a > b ? a : b)

/** 
 * Fortran wrapper for device symmetry apply vector
 */
void opencl_symmetry_apply_vector(void *xmsk, void *ymsk, void *zmsk,
				  void *x, void *y, void *z,
				  int *m, int *n, int *l) {

  cl_int err;
   
  if (symmetry_program == NULL)
    opencl_kernel_jit(symmetry_kernel, (cl_program *) &symmetry_program);
  
  cl_kernel kernel = clCreateKernel(symmetry_program,
				    "symmetry_apply_vector_kernel", &err);
  CL_CHECK(err);
 
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &xmsk));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &ymsk));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &zmsk));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &y));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &z));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(int), m));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), n));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(int), l));
  
  const int max_len = MAX(MAX(*m, *n), *l);
  const size_t global_item_size = 256 * max_len;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				  NULL, &global_item_size, &local_item_size,
				  0, NULL, NULL));
}
