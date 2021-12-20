#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>

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

  switch(*lx) {
  case 2:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx2", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 3:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx3", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 4:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx4", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 5:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx5", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 6:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx6", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 7:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx7", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 8:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx8", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 9:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx9", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 10:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx10", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 11:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx11", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
    break;
  case 12:
    {
      cl_kernel kernel = clCreateKernel(ax_helm_program,
					"ax_helm_kernel_lx12", &err);

      err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w);
      err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u);
      err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx);
      err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy);
      err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz);
      err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dxt);
      err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dyt);
      err = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dzt);
      err = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &h1);
      err = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g11);
      err = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g22);      
      err = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g33);
      err = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g12);
      err = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g13);
      err = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g23);
	
      err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
				   NULL, &global_item_size, &local_item_size,
				   0, NULL, NULL);
    }
  }
}
