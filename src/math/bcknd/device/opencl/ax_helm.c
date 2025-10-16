/*
 Copyright (c) 2021-2025, The Neko Authors
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>
#include <common/neko_log.h>

#include "ax_helm_kernel.cl.h"

int *autotune_ax_helm = NULL;

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

  const size_t global_item_size = 256 * (*nelv);
  const size_t local_item_size = 256;

  size_t global_kstep[2];
  size_t local_kstep[2];
  local_kstep[0] = (*lx);
  local_kstep[1] = (*lx);
  global_kstep[0] = (*nelv) * (*lx);
  global_kstep[1] = (*lx);

  if (autotune_ax_helm == NULL) {
    autotune_ax_helm = malloc(13 * sizeof(int));
    memset(autotune_ax_helm, 0, 13 * sizeof(int));
  }

#define STR(X) #X
#define CASE_1D(LX, QUEUE, EVENT)                                               \
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
      CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) QUEUE,                 \
                                     kernel, 1, NULL, &global_item_size,        \
                                     &local_item_size, 0, NULL, EVENT));        \
                                                                                \
    }


#define CASE_KSTEP(LX, QUEUE, EVENT)                                            \
    {                                                                           \
      cl_kernel kernel = clCreateKernel(ax_helm_program,                        \
                                        STR(ax_helm_kernel_kstep_lx##LX), &err);\
      CL_CHECK(err);                                                            \
                                                                                \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &w));         \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &u));         \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dx));        \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dy));        \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dz));        \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &h1));        \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &g11));       \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &g22));       \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &g33));       \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &g12));       \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g13));      \
      CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g23));      \
                                                                                \
     CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) QUEUE,                  \
                                     kernel, 2, NULL, global_kstep,             \
                                     local_kstep, 0, NULL, EVENT));             \
                                                                                \
    }


#define CASE(LX)                                                                \
  case LX:                                                                      \
      if(autotune_ax_helm[LX] == 0 ) {                                          \
        char *env_value = NULL;                                                 \
        char neko_log_buf[80];                                                  \
        env_value = getenv("NEKO_AUTOTUNE");                                    \
                                                                                \
        sprintf(neko_log_buf, "Autotune Ax helm (lx: %d)", *lx);                \
        log_section(neko_log_buf);                                              \
        if(env_value) {                                                         \
          if( !strcmp(env_value,"1D") ) {                                       \
            CASE_1D(LX, glb_cmd_queue, NULL);                                   \
            sprintf(neko_log_buf,"Set by env : 1 (1D)");                        \
            log_message(neko_log_buf);                                          \
            autotune_ax_helm[LX] = 1;                                           \
          } else if( !strcmp(env_value,"KSTEP") ) {                             \
            CASE_KSTEP(LX, glb_cmd_queue, NULL);                                \
            sprintf(neko_log_buf,"Set by env : 2 (KSTEP)");                     \
            log_message(neko_log_buf);                                          \
            autotune_ax_helm[LX] = 2;                                           \
          } else {                                                              \
            sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");       \
            log_error(neko_log_buf);                                            \
          }                                                                     \
        }                                                                       \
        else {                                                                  \
          CL_CHECK(clFinish(glb_cmd_queue));                                    \
          cl_event perf_event, sync_event;                                      \
          cl_ulong start, end;                                                  \
          CL_CHECK(clEnqueueMarker(glb_cmd_queue, &sync_event));                \
          CL_CHECK(clEnqueueBarrier(prf_cmd_queue));                            \
          CL_CHECK(clEnqueueWaitForEvents(prf_cmd_queue, 1, &sync_event));      \
                                                                                \
          double elapsed1 = 0.0;                                                \
          for(int i = 0; i < 100; i++) {                                        \
            CASE_1D(LX, prf_cmd_queue, &perf_event);                            \
            CL_CHECK(clWaitForEvents(1, &perf_event));                          \
            CL_CHECK(clGetEventProfilingInfo(perf_event,                        \
                                             CL_PROFILING_COMMAND_START,        \
                                             sizeof(cl_ulong), &start, NULL));  \
            CL_CHECK(clGetEventProfilingInfo(perf_event,                        \
                                             CL_PROFILING_COMMAND_END,          \
                                             sizeof(cl_ulong), &end, NULL));    \
            elapsed1 += (end - start)*1.0e-6;                                   \
          }                                                                     \
                                                                                \
          double elapsed2 = 0.0;                                                \
          for(int i = 0; i < 100; i++) {                                        \
            CASE_KSTEP(LX, prf_cmd_queue, &perf_event);                         \
            CL_CHECK(clWaitForEvents(1, &perf_event));                          \
            CL_CHECK(clGetEventProfilingInfo(perf_event,                        \
                                             CL_PROFILING_COMMAND_START,        \
                                             sizeof(cl_ulong), &start, NULL));  \
            CL_CHECK(clGetEventProfilingInfo(perf_event,                        \
                                             CL_PROFILING_COMMAND_END,          \
                                             sizeof(cl_ulong), &end, NULL));    \
            elapsed2 += (end - start)*1.0e-6;                                   \
          }                                                                     \
                                                                                \
          CL_CHECK(clFinish(prf_cmd_queue));                                    \
          CL_CHECK(clEnqueueMarker(prf_cmd_queue, &sync_event));                \
          int krnl_strtgy = (elapsed1 < elapsed2 ? 1 : 2);                      \
          sprintf(neko_log_buf, "Chose      : %d (%s)", krnl_strtgy,            \
                  (krnl_strtgy > 1 ? "KSTEP" : "1D"));                          \
          autotune_ax_helm[LX] = krnl_strtgy;                                   \
          log_message(neko_log_buf);                                            \
          clEnqueueBarrier(glb_cmd_queue);                                      \
          clEnqueueWaitForEvents(glb_cmd_queue, 1, &sync_event) ;               \
        }                                                                       \
        log_end_section();                                                      \
      } else if (autotune_ax_helm[LX] == 1 ) {                                  \
        CASE_1D(LX, glb_cmd_queue, NULL);                                       \
      } else if (autotune_ax_helm[LX] == 2 ) {                                  \
        CASE_KSTEP(LX, glb_cmd_queue, NULL);                                    \
      }                                                                         \
      break

#define CASE_LARGE(LX)                                                          \
    case LX:                                                                    \
      CASE_KSTEP(LX, glb_cmd_queue, NULL);                                      \
      break

  if((*lx) < 13) {
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
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }
  }
  else {
    switch(*lx) {
      CASE_LARGE(13);
      CASE_LARGE(14);
      CASE_LARGE(15);
      CASE_LARGE(16);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }
  }
}

/**
 * Fortran wrapper for device OpenCL Ax (vector version)
 */
void opencl_ax_helm_vector(void *au, void *av, void *aw,
                           void *u, void *v, void *w,
                           void *dx, void *dy, void *dz,
                           void *dxt, void *dyt, void *dzt, void *h1,
                           void *g11, void *g22, void *g33, void *g12,
                           void *g13, void *g23, int *nelv, int *lx) {

  cl_int err;

  if (ax_helm_program == NULL)
      opencl_kernel_jit(ax_helm_kernel, (cl_program *) &ax_helm_program);

  size_t * global_kstep = (size_t *) malloc(sizeof(size_t) * 2);
  size_t * local_kstep = (size_t *) malloc(sizeof(size_t) * 2);
  local_kstep[0] = (*lx);
  local_kstep[1] = (*lx);
  global_kstep[0] = (*nelv) * (*lx);
  global_kstep[1] = (*lx);


#define STR(X) #X
#define CASE_VECTOR(LX)                                                         \
  case LX:                                                                      \
    {                                                                           \
      cl_kernel kernel =                                                        \
        clCreateKernel(ax_helm_program,                                         \
                       STR(ax_helm_kernel_vector_kstep_lx##LX), &err);          \
      CL_CHECK(err);                                                            \
                                                                                \
      CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &au));        \
      CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &av));        \
      CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &aw));        \
      CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &u));         \
      CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &v));         \
      CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &w));         \
      CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dx));        \
      CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dy));        \
      CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &dz));        \
      CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &h1));        \
      CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &g11));      \
      CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &g22));      \
      CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &g33));      \
      CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &g12));      \
      CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &g13));      \
      CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &g23));      \
                                                                                \
     CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue,          \
                                     kernel, 2, NULL, global_kstep,             \
                                     local_kstep, 0, NULL, NULL));              \
                                                                                \
    }                                                                           \
    break

  switch(*lx) {
    CASE_VECTOR(2);
    CASE_VECTOR(3);
    CASE_VECTOR(4);
    CASE_VECTOR(5);
    CASE_VECTOR(6);
    CASE_VECTOR(7);
    CASE_VECTOR(8);
    CASE_VECTOR(9);
    CASE_VECTOR(10);
    CASE_VECTOR(11);
    CASE_VECTOR(12);
    CASE_VECTOR(13);
    CASE_VECTOR(14);
    CASE_VECTOR(15);
    CASE_VECTOR(16);
  default:
    {
      fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
      exit(1);
    }
  }

  free(global_kstep);
  free(local_kstep);
}
