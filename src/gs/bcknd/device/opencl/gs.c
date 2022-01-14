#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <limits.h>
#include <stdio.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>

#include "gs_kernels.cl.h"

#define GS_OP_ADD  1
#define GS_OP_MUL  2
#define GS_OP_MIN  3
#define GS_OP_MAX  4

/** 
 * Fortran wrapper for device gather kernels
 */
void opencl_gather_kernel(void *v, int *m, int *o, void *dg,
                          void *u, int *n, void *gd, int *nb,
                          void *b, void *bo, int *op) {
  cl_int err;
  
  if (gs_program == NULL)
    opencl_kernel_jit(gs_kernels, (cl_program *) &gs_program);
  
  const int nblks = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nblks;
  const size_t local_item_size = 256;
  
  
  switch (*op) {
  case GS_OP_ADD:
    {
      const real zero = 0;
      CL_CHECK(clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
                                   v, &zero, sizeof(real), 0,
                                   (*m) * sizeof(real), 0, NULL, NULL));
            
      cl_kernel kernel = clCreateKernel(gs_program,
                                        "gather_kernel_add", &err);
      CL_CHECK(err);
  
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v));
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(int), m));
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), o));
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dg));
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &u));
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &gd));
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), nb));
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &b));
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bo));
    
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                      1, NULL, &global_item_size,
                                      &local_item_size, 0, NULL, NULL));
    }
    break;
  case GS_OP_MUL:
    {
      const real one = 0;
      CL_CHECK(clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
                                   v, &one, sizeof(real), 0,
                                   (*m) * sizeof(real), 0, NULL, NULL));
      
      cl_kernel kernel = clCreateKernel(gs_program,
                                        "gather_kernel_mul", &err);
      CL_CHECK(err);
      
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v));
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(int), m));
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), o));
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dg));
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &u));
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &gd));
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), nb));
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &b));
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bo));
    
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                      1, NULL, &global_item_size,
                                      &local_item_size, 0, NULL, NULL));
    }
    break;
  case GS_OP_MIN:
    {
      const real rmax = (real) INT_MAX;
      
      CL_CHECK(clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
                                   v, &rmax, sizeof(real), 0,
                                   (*m) * sizeof(real), 0, NULL, NULL));
      
      cl_kernel kernel = clCreateKernel(gs_program,
                                        "gather_kernel_min", &err);
      CL_CHECK(err);
        
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v));
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(int), m));
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), o));
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dg));
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &u));
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &gd));
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), nb));
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &b));
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bo));
    
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                      1, NULL, &global_item_size,
                                      &local_item_size, 0, NULL, NULL));
    }
    break;
  case GS_OP_MAX:
    {
      const real rmin = (real) -INT_MAX;
      
      CL_CHECK(clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
                                   v, &rmin, sizeof(real), 0,
                                   (*m) * sizeof(real), 0, NULL, NULL));
      
      cl_kernel kernel = clCreateKernel(gs_program,
                                        "gather_kernel_max", &err);
      CL_CHECK(err);
      
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v));
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(int), m));
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(int), o));
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dg));
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &u));
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), n));
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &gd));
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(int), nb));
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &b));
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bo));
    
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel,
                                      1, NULL, &global_item_size,
                                      &local_item_size, 0, NULL, NULL));
    }
    break;
  }
}

/**
 * Fortran wrapper for device scatter kernel
 */
void opencl_scatter_kernel(void *v, int *m, void *dg,
                           void *u, int *n, void *gd,
                           int *nb, void *b, void *bo) {
  cl_int err;

  if (gs_program == NULL)
    opencl_kernel_jit(gs_kernels, (cl_program *) &gs_program);
  
  cl_kernel kernel = clCreateKernel(gs_program, "scatter_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(int), m));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dg));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &u));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(int), n));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &gd));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(int), nb));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &b));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &bo));
  
  const int nblks = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nblks;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));
}
