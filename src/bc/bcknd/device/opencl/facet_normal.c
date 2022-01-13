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

#include "facet_normal_kernel.cl.h"

/** 
 * Fortran wrapper for device facet normal apply surfvec
 */
void opencl_facet_normal_apply_surfvec(void *msk, void *facet,
                                       void *x, void *y, void *z,
                                       void *u, void *v, void *w,
                                       void *nx, void * ny, void *nz,
                                       void *area, int *lx, int *m) {

  cl_int err;
  
  if (facet_normal_program == NULL)
    opencl_kernel_jit(facet_normal_kernel, (cl_program *) &facet_normal_program);
  
  cl_kernel kernel = clCreateKernel(facet_normal_program,
                                    "facet_normal_apply_surfvec_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &msk));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &facet));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &x));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &y));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &z));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &w));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &nx));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &ny));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &nz));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &facet));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(int), lx));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(int), m));
  
  const int nb = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}
