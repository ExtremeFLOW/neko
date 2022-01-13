#ifndef __CL_CHECK_H
#define __CL_CHECK_H

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

void opencl_check(const char *fname, const int line, const cl_int err);

#define CL_CHECK(err) opencl_check(__FILE__, __LINE__, err)

#endif
