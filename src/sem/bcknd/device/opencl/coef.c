/*
 Copyright (c) 2022-2023, The Neko Authors
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

#include <stdio.h>
#include <math.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>

#include "coef_kernel.cl.h"

/** 
 * Fortran wrapper for generating geometric factors
 */
void opencl_coef_generate_geo(void *G11, void *G12, void *G13, 
                              void *G22, void *G23, void *G33, 
                              void *drdx, void *drdy, void *drdz,
                              void *dsdx, void *dsdy, void *dsdz, 
                              void *dtdx, void *dtdy, void *dtdz, 
                              void *jacinv, void *w3, int *nel, 
                              int *lx, int *gdim) {

  cl_int err;
  int i;
  if (coef_program == NULL)
    opencl_kernel_jit(coef_kernel, (cl_program *) &coef_program);
  
  const size_t global_item_size = 256 * (*nel);
  const size_t local_item_size = 256;

#define STR(X) #X
#define GEO_CASE(LX)                                                            \
  case LX:                                                                      \
    {                                                                           \
      cl_kernel kernel = clCreateKernel(coef_program,                           \
                                   STR(coef_generate_geo_kernel_lx##LX), &err); \
      CL_CHECK(err);                                                            \
                                                                                \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &G11));       \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &G12));       \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &G13));       \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &G22));       \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &G23));       \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &G33));       \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &drdx));      \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &drdy));      \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &drdz));      \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &dsdx));      \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &dsdy));     \
      CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &dsdz));     \
      CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &dtdx));     \
      CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &dtdy));     \
      CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &dtdz));     \
      CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &jacinv));   \
      CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &w3));       \
      CL_CHECK(clSetKernelArg(kernel, 17, sizeof(int), gdim));                  \
                                                                                \
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,         \
                                      kernel, 1, NULL, &global_item_size,       \
                                      &local_item_size, 0, NULL, NULL));        \
    }                                                                           \
    break
    
  switch(*lx) {
    GEO_CASE(2);
    GEO_CASE(3);
    GEO_CASE(4);
    GEO_CASE(5);
    GEO_CASE(6);
    GEO_CASE(7);
    GEO_CASE(8);
    GEO_CASE(9);
    GEO_CASE(10);
    GEO_CASE(11);
    GEO_CASE(12);
    GEO_CASE(13);
    GEO_CASE(14);
    GEO_CASE(15);
    GEO_CASE(16);
  }
}

/** 
 * Fortran wrapper for generating geometric factors
 */
void opencl_coef_generate_dxyzdrst(void *drdx, void *drdy, void *drdz, 
                                   void *dsdx, void *dsdy, void *dsdz, 
                                   void *dtdx, void *dtdy, void *dtdz, 
                                   void *dxdr, void *dydr, void *dzdr, 
                                   void *dxds, void *dyds, void *dzds, 
                                   void *dxdt, void *dydt, void *dzdt,
                                   void *dx, void *dy, void *dz, 
                                   void *x, void *y, void *z,
                                   void *jacinv, void *jac,
                                   int *lx, int *nel)  {

  cl_int err;
  int i;
  if (coef_program == NULL)
    opencl_kernel_jit(coef_kernel, (cl_program *) &coef_program);
  
  const int n = (*nel) * (*lx) * (*lx) * (*lx);
  const size_t global_item_size_dxyz = 256 * (*nel);
  const size_t global_item_size_drst = 256 * n;
  const size_t local_item_size = 256;

#define STR(X) #X
#define DXYZDRST_CASE(LX)                                                       \
  case LX:                                                                      \
    {                                                                           \
      cl_kernel kernel = clCreateKernel(coef_program,                           \
                                  STR(coef_generate_dxyz_kernel_lx##LX), &err); \
      CL_CHECK(err);                                                            \
                                                                                \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &dxdr));      \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &dydr));      \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dzdr));      \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dxds));      \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dyds));      \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dzds));      \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dxdt));      \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dydt));      \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &dzdt));      \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &dx));        \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &dy));       \
      CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &dz));       \
      CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &x));        \
      CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &y));        \
      CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &z));        \
                                                                                \
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,         \
                                      kernel, 1, NULL, &global_item_size_dxyz,  \
                                      &local_item_size, 0, NULL, NULL));        \
    }                                                                           \
    break
    
  switch(*lx) {
    DXYZDRST_CASE(2);
    DXYZDRST_CASE(3);
    DXYZDRST_CASE(4);
    DXYZDRST_CASE(5);
    DXYZDRST_CASE(6);
    DXYZDRST_CASE(7);
    DXYZDRST_CASE(8);
    DXYZDRST_CASE(9);
    DXYZDRST_CASE(10);
    DXYZDRST_CASE(11);
    DXYZDRST_CASE(12);
    DXYZDRST_CASE(13);
    DXYZDRST_CASE(14);
    DXYZDRST_CASE(15);
    DXYZDRST_CASE(16);
  }

  cl_kernel kernel = clCreateKernel(coef_program,
                                    "coef_generate_drst_kernel", &err);
  CL_CHECK(err);                                                            
                                                                                
  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &jac));      
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &jacinv));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &drdx));      
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &drdy));      
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &drdz));      
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dsdx));      
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dsdy));      
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dsdz));      
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &dtdx));      
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &dtdy));      
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &dtdz));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &dxdr));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &dydr));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &dzdr));
  CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &dxds));
  CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &dyds));
  CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), (void *) &dzds));
  CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), (void *) &dxdt));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), (void *) &dydt));
  CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), (void *) &dzdt));     
  CL_CHECK(clSetKernelArg(kernel, 20, sizeof(int), &n));        
  
  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,         
                                  kernel, 1, NULL, &global_item_size_drst,  
                                  &local_item_size, 0, NULL, NULL));   
}
