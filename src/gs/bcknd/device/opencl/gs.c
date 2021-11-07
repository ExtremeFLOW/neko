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
      err = clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
				v, &zero, sizeof(real), 0,
				(*m) * sizeof(real), 0, NULL, NULL);
	    
      cl_kernel kernel = clCreateKernel(gs_program,
					"gather_kernel_add", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v);
      err = clSetKernelArg(kernel, 1, sizeof(int), m);
      err = clSetKernelArg(kernel, 2, sizeof(int), o);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dg);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 5, sizeof(int), n);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &gd);
      err = clSetKernelArg(kernel, 7, sizeof(int), nb);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &b);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bo);
    
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case GS_OP_MUL:
    {
      const real one = 0;
      err = clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
				v, &one, sizeof(real), 0,
				(*m) * sizeof(real), 0, NULL, NULL);
      
      cl_kernel kernel = clCreateKernel(gs_program,
					"gather_kernel_mul", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v);
      err = clSetKernelArg(kernel, 1, sizeof(int), m);
      err = clSetKernelArg(kernel, 2, sizeof(int), o);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dg);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 5, sizeof(int), n);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &gd);
      err = clSetKernelArg(kernel, 7, sizeof(int), nb);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &b);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bo);
    
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case GS_OP_MIN:
    {
      const real rmax = (real) INT_MAX;
      
      err = clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
				v, &rmax, sizeof(real), 0,
				(*m) * sizeof(real), 0, NULL, NULL);

      cl_kernel kernel = clCreateKernel(gs_program,
					"gather_kernel_min", &err);
      
      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v);
      err = clSetKernelArg(kernel, 1, sizeof(int), m);
      err = clSetKernelArg(kernel, 2, sizeof(int), o);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dg);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 5, sizeof(int), n);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &gd);
      err = clSetKernelArg(kernel, 7, sizeof(int), nb);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &b);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bo);
    
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case GS_OP_MAX:
    {
      const real rmin = (real) -INT_MAX;
      
      err = clEnqueueFillBuffer((cl_command_queue) glb_cmd_queue,
				v, &rmin, sizeof(real), 0,
				(*m) * sizeof(real), 0, NULL, NULL);
      
      cl_kernel kernel = clCreateKernel(gs_program,
					"gather_kernel_max", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v);
      err = clSetKernelArg(kernel, 1, sizeof(int), m);
      err = clSetKernelArg(kernel, 2, sizeof(int), o);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dg);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 5, sizeof(int), n);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &gd);
      err = clSetKernelArg(kernel, 7, sizeof(int), nb);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &b);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &bo);
    
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
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

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &v);
  err = clSetKernelArg(kernel, 1, sizeof(int), m);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dg);
  err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &u);
  err = clSetKernelArg(kernel, 4, sizeof(int), n);
  err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &gd);
  err = clSetKernelArg(kernel, 6, sizeof(int), nb);
  err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &b);
  err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &bo);
  
  const int nblks = ((*m) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nblks;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);
}
