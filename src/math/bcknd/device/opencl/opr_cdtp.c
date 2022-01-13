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

#include "cdtp_kernel.cl.h"

/** 
 * Fortran wrapper for device OpenCL \f$ D^T X \f$
 */
void opencl_cdtp(void *dtx, void *x,
                 void *dr, void *ds, void *dt,
                 void *dxt, void *dyt, void *dzt,
                 void *B, void *jac, int *nel, int *lx) {
  cl_int err;
  
  if (cdtp_program == NULL)
    opencl_kernel_jit(cdtp_kernel, (cl_program *) &cdtp_program);
  
  const int nb = ((*nel) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

#define STR(X) #X
#define CASE(LX)                                                                \
  case LX:                                                                      \
    {                                                                           \
      cl_kernel kernel = clCreateKernel(cdtp_program,                           \
                                        STR(cdtp_kernel_lx##LX), &err);         \
      CL_CHECK(err)                                                             \
                                                                                \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &dtx))        \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &x))          \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr))         \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds))         \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt))         \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt))        \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt))        \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt))        \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &B))          \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &jac))        \
                                                                                \
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,         \
                                      kernel, 1, NULL, &global_item_size,       \
                                      &local_item_size, 0, NULL, NULL))         \
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

