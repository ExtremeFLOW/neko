#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <device/device_config.h>
#include <device/opencl/check.h>


/**
 * JIT compile OpenCL kernels
 */
void opencl_kernel_jit(const char *kernel, cl_program *program)
{
  cl_int err;

  (*program) = clCreateProgramWithSource((cl_context) glb_ctx, 1,
                                           &kernel, NULL, &err);
  CL_CHECK(err);
  
  err = clBuildProgram((*program), 1, (cl_device_id *) &glb_device_id,
                       NULL, NULL, NULL);

  if (err == CL_BUILD_PROGRAM_FAILURE) {
    size_t log_size;
    clGetProgramBuildInfo((*program), glb_device_id,
                          CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    char *log = (char *) malloc(log_size);

    clGetProgramBuildInfo((*program), glb_device_id,
                          CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

    fprintf(stderr, "%s\n", log);

  }
  
}


