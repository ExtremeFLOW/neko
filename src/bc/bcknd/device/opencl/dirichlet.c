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

#include "dirichlet_kernel.cl.h"

/** 
 * Fortran wrapper for device dirichlet apply scalar
 */
void opencl_dirichlet_apply_scalar(void *msk, void *x,
                                 real *g, int *m) {
  cl_int err;
  
  if (dirichlet_program == NULL)
    opencl_kernel_jit(dirichlet_kernel, (cl_program *) &dirichlet_program);
  
  cl_kernel kernel = clCreateKernel(dirichlet_program,
                                    "dirichlet_apply_scalar_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &msk));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(real), g));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(int), m));
  
  const int nb = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
}

/** 
 * Fortran wrapper for device dirichlet apply vector
 */
void opencl_dirichlet_apply_vector(void *msk, void *x, void *y,
                                 void *z, real *g, int *m) {
  cl_int err;
  
  if (dirichlet_program == NULL)
    opencl_kernel_jit(dirichlet_kernel, (cl_program *) &dirichlet_program);
  
  cl_kernel kernel = clCreateKernel(dirichlet_program,
                                    "dirichlet_apply_vector_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &msk));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &y));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &z));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(real), g));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), m));
  
  const int nb = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}
