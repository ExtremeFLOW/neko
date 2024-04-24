/*
 Copyright (c) 2022-2024, The Neko Authors
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
#include <stdarg.h>
#include <stdlib.h>
#include <device/opencl/check.h>

#define CL_ERR_STR(err) case err: return #err

/**
 * Convert an OpenCL error to a string
 */
const char *clGetError(const cl_int err) {
  
  switch(err) {
    CL_ERR_STR(CL_SUCCESS);
    CL_ERR_STR(CL_DEVICE_NOT_FOUND);
    CL_ERR_STR(CL_DEVICE_NOT_AVAILABLE);
    CL_ERR_STR(CL_COMPILER_NOT_AVAILABLE);
    CL_ERR_STR(CL_MEM_OBJECT_ALLOCATION_FAILURE);
    CL_ERR_STR(CL_OUT_OF_RESOURCES);
    CL_ERR_STR(CL_OUT_OF_HOST_MEMORY);
    CL_ERR_STR(CL_PROFILING_INFO_NOT_AVAILABLE);
    CL_ERR_STR(CL_MEM_COPY_OVERLAP);
    CL_ERR_STR(CL_IMAGE_FORMAT_MISMATCH);
    CL_ERR_STR(CL_IMAGE_FORMAT_NOT_SUPPORTED);
    CL_ERR_STR(CL_BUILD_PROGRAM_FAILURE);
    CL_ERR_STR(CL_MAP_FAILURE);
    CL_ERR_STR(CL_MISALIGNED_SUB_BUFFER_OFFSET);
    CL_ERR_STR(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST);
    CL_ERR_STR(CL_COMPILE_PROGRAM_FAILURE);
    CL_ERR_STR(CL_LINKER_NOT_AVAILABLE);
    CL_ERR_STR(CL_LINK_PROGRAM_FAILURE);
    CL_ERR_STR(CL_DEVICE_PARTITION_FAILED);
    CL_ERR_STR(CL_KERNEL_ARG_INFO_NOT_AVAILABLE);
    CL_ERR_STR(CL_INVALID_VALUE);
    CL_ERR_STR(CL_INVALID_DEVICE_TYPE);
    CL_ERR_STR(CL_INVALID_PLATFORM);
    CL_ERR_STR(CL_INVALID_DEVICE);
    CL_ERR_STR(CL_INVALID_CONTEXT);
    CL_ERR_STR(CL_INVALID_QUEUE_PROPERTIES);
    CL_ERR_STR(CL_INVALID_COMMAND_QUEUE);               
    CL_ERR_STR(CL_INVALID_HOST_PTR);
    CL_ERR_STR(CL_INVALID_MEM_OBJECT);                       
    CL_ERR_STR(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR);          
    CL_ERR_STR(CL_INVALID_IMAGE_SIZE);                       
    CL_ERR_STR(CL_INVALID_SAMPLER);                          
    CL_ERR_STR(CL_INVALID_BINARY);                           
    CL_ERR_STR(CL_INVALID_BUILD_OPTIONS);                    
    CL_ERR_STR(CL_INVALID_PROGRAM);                
    CL_ERR_STR(CL_INVALID_PROGRAM_EXECUTABLE);               
    CL_ERR_STR(CL_INVALID_KERNEL_NAME);                      
    CL_ERR_STR(CL_INVALID_KERNEL_DEFINITION);                
    CL_ERR_STR(CL_INVALID_KERNEL);                           
    CL_ERR_STR(CL_INVALID_ARG_INDEX);                        
    CL_ERR_STR(CL_INVALID_ARG_VALUE);                        
    CL_ERR_STR(CL_INVALID_ARG_SIZE);                         
    CL_ERR_STR(CL_INVALID_KERNEL_ARGS);                      
    CL_ERR_STR(CL_INVALID_WORK_DIMENSION);                   
    CL_ERR_STR(CL_INVALID_WORK_GROUP_SIZE);                  
    CL_ERR_STR(CL_INVALID_WORK_ITEM_SIZE);
    CL_ERR_STR(CL_INVALID_GLOBAL_OFFSET);                    
    CL_ERR_STR(CL_INVALID_EVENT_WAIT_LIST);                  
    CL_ERR_STR(CL_INVALID_EVENT);                            
    CL_ERR_STR(CL_INVALID_OPERATION);                        
    CL_ERR_STR(CL_INVALID_GL_OBJECT);                        
    CL_ERR_STR(CL_INVALID_BUFFER_SIZE);                      
    CL_ERR_STR(CL_INVALID_MIP_LEVEL);
    CL_ERR_STR(CL_INVALID_GLOBAL_WORK_SIZE);
    CL_ERR_STR(CL_INVALID_PROPERTY);
    CL_ERR_STR(CL_INVALID_IMAGE_DESCRIPTOR);
    CL_ERR_STR(CL_INVALID_COMPILER_OPTIONS);
    CL_ERR_STR(CL_INVALID_LINKER_OPTIONS);
    CL_ERR_STR(CL_INVALID_DEVICE_PARTITION_COUNT);
  }
}
  

/**
 * Check a OpenCL return code
 */
void opencl_check(const char *fname, const int line, const cl_int err)
{
  if (err != CL_SUCCESS) {
    const char *err_str = clGetError(err);
    fprintf(stderr, "%s in %s:%d \n", err_str, fname, line);
    exit(1);
  }						  
}
