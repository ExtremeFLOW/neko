/*
 Copyright (c) 2022, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>

#include "gmres_kernel.cl.h"

real *gmres_bf1 = NULL;
cl_mem gmres_bfd1 = NULL;

/** 
 * Fortran wrapper for device gmres part2
 */
real opencl_gmres_part2(void *w, void *v, void *h,
                        void *mult, int *j, int *n) {
  cl_int err;
  cl_event kern_wait;
  
  if (gmres_program == NULL)
    opencl_kernel_jit(gmres_kernel, (cl_program *) &gmres_program);

  const int nb = ((*n) + 256 - 1) / 256; 
  const size_t global_item_size = 256 * nb;                                    
  const size_t local_item_size = 256;

  if (gmres_bf1 == NULL){
    gmres_bf1 = (real *) malloc(nb * sizeof(real));
    gmres_bfd1 = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE,                    
                                nb * sizeof(real), NULL, &err);
  }
      
  cl_kernel kernel = clCreateKernel(gmres_program, "gmres_part2_kernel", &err);
  CL_CHECK(err);                                                           
                                                                               
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w));            
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &v));      
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &h));      
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &gmres_bfd1));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(int), j));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(int), n));      
    
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,        
                                  kernel, 1, NULL, &global_item_size,      
                                  &local_item_size, 0, NULL, &kern_wait));

  CL_CHECK(clEnqueueReadBuffer((cl_command_queue) glb_cmd_queue,             
                               gmres_bfd1, CL_TRUE, 0, nb * sizeof(real),      
                               gmres_bf1, 1, &kern_wait, NULL));

  real res1 = 0.0;
  for (int i = 0; i < nb; i++) {
    res1 += gmres_bf1[i];
  }
  return res1;
}
