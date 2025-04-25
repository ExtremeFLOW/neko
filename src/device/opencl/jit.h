#ifndef __CL_JIT_H
#define __CL_JIT_H

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

void opencl_kernel_jit(const char *kernel, cl_program *program);

#endif
