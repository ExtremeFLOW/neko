/*

  Fortran wrappers for various CUDA calls
  
*/

#include <cstdio>

#define CUDA_MEMCPY_HTOD 1
#define CUDA_MEMCPY_DTOH 2

#define CHECK_ERROR(cuda_call)                        \
  cudaError_t err = (cuda_call);                      \
  if (err != cudaSuccess) {			      \
    fprintf(stderr, "%s\n", cudaGetErrorString(err)); \
    return 1;                                         \
  }                                                   \


extern "C" {

  int cudaMalloc_wrapper(void **ptr, size_t *size) {
    CHECK_ERROR(cudaMalloc(ptr, *size));
    return 0;   
  }

  int cudaFree_wrapper(void *ptr) {
    CHECK_ERROR(cudaFree(ptr));
    return 0;
  }

  int cudaMemcpy_wrapper(void *ptr_h, void *ptr_d, size_t *size, int *dir) {
    if ((*dir) == CUDA_MEMCPY_HTOD) {
      CHECK_ERROR(cudaMemcpy(ptr_d, ptr_h, *size, cudaMemcpyHostToDevice));
    }
    else if ((*dir) == CUDA_MEMCPY_DTOH) {
      CHECK_ERROR(cudaMemcpy(ptr_h, ptr_d, *size, cudaMemcpyDeviceToHost));
    }
    else {
      fprintf(stderr, "Invalid direction\n");
      return 1;
    }
    return 0;
  }

  int cudaDeviceSynchronize_wrapper() {
    CHECK_ERROR(cudaDeviceSynchronize());
    return 0;
  }
  
}
