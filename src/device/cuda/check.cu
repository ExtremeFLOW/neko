#include <stdio.h>
#include <device/cuda/check.h>

/**
 * Check a CUDA return code
 */
void cuda_check(const char *fname, const int line, const cudaError_t err)
{
  if (err != cudaSuccess) {
    fprintf(stderr, "%s in %s:%d \n", cudaGetErrorString(err), fname, line);
  }						  
}
