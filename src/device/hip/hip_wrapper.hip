/*

  Fortran wrappers for various HIP calls
  
*/

#include <hip/hip_runtime.h>
#include <cstdio>

#define HIP_MEMCPY_HTOD 1
#define HIP_MEMCPY_DTOH 2

#define CHECK_ERROR(hip_call)                        \
  hipError_t err = (hip_call);                       \
  if (err != hipSuccess) {			     \
    fprintf(stderr, "%s\n", hipGetErrorString(err)); \
    return 1;                                        \
  }                                                  \


extern "C" {

  int hipMalloc_wrapper(void **ptr, size_t *size) {
    CHECK_ERROR(hipMalloc(ptr, *size));
    return 0;   
  }

  int hipFree_wrapper(void *ptr) {
    CHECK_ERROR(hipFree(ptr));
    return 0;
  }

  int hipMemcpy_wrapper(void *ptr_h, void *ptr_d, size_t *size, int *dir) {
    if ((*dir) == HIP_MEMCPY_HTOD) {
      CHECK_ERROR(hipMemcpy(ptr_d, ptr_h, *size, hipMemcpyHostToDevice));
    }
    else if ((*dir) == HIP_MEMCPY_DTOH) {
      CHECK_ERROR(hipMemcpy(ptr_h, ptr_d, *size, hipMemcpyDeviceToHost));
    }
    else {
      fprintf(stderr, "Invalid direction\n");
      return 1;
    }
    return 0;
  }

  int hipDeviceSynchronize_wrapper() {
    CHECK_ERROR(hipDeviceSynchronize());
    return 0;
  }
  
}
