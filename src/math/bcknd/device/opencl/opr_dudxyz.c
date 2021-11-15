#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>

#include "dudxyz_kernel.cl.h"

/** 
 * Fortran wrapper for device OpenCL derivative kernels
 */
void opencl_dudxyz(void *du, void *u,
		   void *dr, void *ds, void *dt,
		   void *dx, void *dy, void *dz,
		   void *jacinv, int *nel, int *lx) {
  cl_int err;
  
  if (dudxyz_program == NULL)
    opencl_kernel_jit(dudxyz_kernel, (cl_program *) &dudxyz_program);
  
  const int nb = ((*nel) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;
  
  switch(*lx) {
  case 2:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx2", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;     
  case 3:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx3", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;     
  case 4:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx4", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;     
  case 5:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx5", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;     
  case 6:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx6", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;     
  case 7:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx7", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;     
  case 8:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx8", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;     
  case 9:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx9", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;     
  case 10:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx10", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }      
    break;
  case 11:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx11", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }      
    break;
  case 12:
    {
      cl_kernel kernel = clCreateKernel(dudxyz_program,
					"dudxyz_kernel_lx12", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &du);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dr);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &ds);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dt);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &jacinv);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }      
    break;
  }
} 

