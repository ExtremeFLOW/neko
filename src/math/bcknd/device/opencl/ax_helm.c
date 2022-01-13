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

#include "ax_helm_kernel.cl.h"

/** 
 * Fortran wrapper for device OpenCL Ax
 */
void opencl_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                    void *dxt, void *dyt, void *dzt, void *h1,
                    void *g11, void *g22, void *g33, void *g12,
                    void *g13, void *g23, int *nelv, int *lx) {

  cl_int err;
  
  if (ax_helm_program == NULL)
      opencl_kernel_jit(ax_helm_kernel, (cl_program *) &ax_helm_program);
  
  const int nb = ((*nelv) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;    

#define STR(X) #X
#define CASE(LX)                                                                \
  case LX:                                                                      \
    {                                                                           \
      cl_kernel kernel = clCreateKernel(ax_helm_program,                        \
                                        STR(ax_helm_kernel_lx##LX), &err);      \
      CL_CHECK(err);                                                            \
                                                                                \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w));         \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u));         \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx));        \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy));        \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz));        \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt));       \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt));       \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt));       \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1));        \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11));       \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22));      \
      CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33));      \
      CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12));      \
      CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13));      \
      CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23));      \
                                                                                \
     CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,          \
                                     kernel, 1, NULL, &global_item_size,        \
                                     &local_item_size, 0, NULL, NULL));         \
                                                                                \
    }                                                                           \
    break

  switch(*lx) {
    CASE(2);
    CASE(3);
    CASE(4);
    CASE(5);
    CASE(6);
    CASE(7);
    CASE(8);
    CASE(9);
    CASE(10);
    CASE(11);
    CASE(12);
  }
}
