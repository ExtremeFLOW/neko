#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>

#include "mathops_kernel.cl.h"

/** Fortran wrapper for opchsign \f$ a = -a \f$ */
void opencl_opchsign(void *a1, void *a2, void *a3, int *gdim, int *n) {    
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opchsign_kernel", &err);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3);
  err = clSetKernelArg(kernel, 3, sizeof(int), gdim);
  err = clSetKernelArg(kernel, 4, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);
}

/** Fortran wrapper for opcolv \f$ a = a * c \f$ */
void opencl_opcolv(void *a1, void *a2, void *a3, void *c, int *gdim, int *n) {
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opcolv_kernel", &err);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3);
  err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &c);
  err = clSetKernelArg(kernel, 4, sizeof(int), gdim);
  err = clSetKernelArg(kernel, 5, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);   
}

/** Fortran wrapper for opcolv3c \f$ a(i) = b(i) * c(i) * d \f$ */
void opencl_opcolv3c(void *a1, void *a2, void *a3,
		     void *b1, void *b2, void *b3,
		     void *c, real *d, int *gdim, int *n) {
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opcolv3c_kernel", &err);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3);
  err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &b1);
  err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &b2);
  err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &b3);
  err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &c);
  err = clSetKernelArg(kernel, 7, sizeof(real), d);
  err = clSetKernelArg(kernel, 8, sizeof(int), gdim);
  err = clSetKernelArg(kernel, 9, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);   
}

/** Fortran wrapper for opadd2cm \f$ a(i) = a + b(i) * c \f$ */
void opencl_opadd2cm(void *a1, void *a2, void *a3, 
		     void *b1, void *b2, void *b3,
		     real *c, int *gdim, int *n) {
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opadd2cm_kernel", &err);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3);
  err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &b1);
  err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &b2);
  err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &b3);
  err = clSetKernelArg(kernel, 6, sizeof(real), c);
  err = clSetKernelArg(kernel, 7, sizeof(int), gdim);
  err = clSetKernelArg(kernel, 8, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);  
}

/** Fortran wrapper for opadd2col \f$ a(i) = a + b(i) * c(i) \f$ */
void opencl_opadd2col(void *a1, void *a2, void *a3, 
		      void *b1, void *b2, void *b3,
		      void *c, int *gdim, int *n) {
  cl_int err;

  if (mathops_program == NULL)
    opencl_kernel_jit(mathops_kernel, (cl_program *) &mathops_program);
  
  cl_kernel kernel = clCreateKernel(mathops_program, "opadd2col_kernel", &err);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &a1);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &a2);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &a3);
  err = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &b1);
  err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &b2);
  err = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &b3);
  err = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &c);
  err = clSetKernelArg(kernel, 7, sizeof(int), gdim);
  err = clSetKernelArg(kernel, 8, sizeof(int), n);
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  err = clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
			       NULL, &global_item_size, &local_item_size,
			       0, NULL, NULL);
}
