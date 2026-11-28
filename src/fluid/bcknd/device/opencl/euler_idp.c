/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of the authors nor the names of its contributors may be
     used to endorse or promote products derived from this software without
     specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
*/

#include <pthread.h>
#include <device/opencl/check.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <fluid/bcknd/device/euler_idp_kernels.cl.h>
#include "../euler_idp_args.h"

static void idp_launch(const char *name, int groups, int width, const idp_arg *args,
                       unsigned count) {
  if (!groups)
    return;
  static pthread_mutex_t program_lock = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_lock(&program_lock);
  if (!euler_idp_program)
    opencl_kernel_jit(euler_idp_kernels, (cl_program *)&euler_idp_program);
  pthread_mutex_unlock(&program_lock);
  cl_int err;
  cl_kernel kernel = clCreateKernel(euler_idp_program, name, &err);
  CL_CHECK(err);
  for (unsigned i = 0; i < count; ++i) {
    switch (args[i].type) {
    case IDP_BUFFER:
      CL_CHECK(clSetKernelArg(kernel, i, sizeof(cl_mem), &args[i].buffer));
      break;
    case IDP_INTEGER:
      CL_CHECK(clSetKernelArg(kernel, i, sizeof(int), &args[i].integer));
      break;
    case IDP_SCALAR:
      CL_CHECK(clSetKernelArg(kernel, i, sizeof(real), &args[i].scalar));
      break;
    }
  }
  size_t local = width, global = (size_t)groups * width;
  CL_CHECK(clEnqueueNDRangeKernel(glb_cmd_queue, kernel, 1, NULL, &global, &local, 0,
                                  NULL, NULL));
  CL_CHECK(clReleaseKernel(kernel));
}
