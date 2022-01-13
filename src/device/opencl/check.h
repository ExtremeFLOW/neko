#ifndef __CL_CHECK_H
#define __CL_CHECK_H

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

const char *opencl_check(const cl_int  err);

#define CL_CHECK(err) opencl_check(err);

#endif
