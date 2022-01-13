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

#include "jacobi_kernel.cl.h"

/** 
 * Fortran wrapper for device pc jacobi
 */
void opencl_jacobi_update(void *d,
			  void *dxt, void *dyt, void *dzt,
			  void *G11, void *G22, void *G33,
			  void *G12, void *G13, void *G23,
			  int *nel, int *lx) {
  cl_int err;
   
  if (jacobi_program == NULL)
    opencl_kernel_jit(jacobi_kernel, (cl_program *) &jacobi_program);

  const int nb = (((*nel) * (*lx) * (*lx) * (*lx)) + 256 - 1) / 256; 
  const size_t global_item_size = 256 * nb;				       
  const size_t local_item_size = 256;				       
    
#define STR(X) #X
#define CASE(LX)                                                               \
  case LX:                                                                     \
    {                                                                          \
      cl_kernel kernel = clCreateKernel(jacobi_program,			       \
					STR(jacobi_kernel_lx##LX), &err);      \
      CL_CHECK(err);                                                           \
                                                                               \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &d));	       \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &dxt));      \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dyt));      \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dzt));      \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &G11));      \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &G22));      \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &G33));      \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &G12));      \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &G13));      \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &G23));      \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(int), nel));		       \
                                                                               \
                                                                               \
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,        \
				      kernel, 1, NULL, &global_item_size,      \
				      &local_item_size, 0, NULL, NULL));       \
    }                                                                          \
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
    CASE(13);
    CASE(14);
    CASE(15);
    CASE(16);
  }
}
