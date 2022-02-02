#include <stdio.h>
#include <hip/hip_runtime.h>
#include <device/hip/check.h>

/**
 * Check a HIP return code
 */
void hip_check(const char *fname, const int line, const hipError_t err)
{
  if (err != hipSuccess) {
    fprintf(stderr, "%s in %s:%d \n", hipGetErrorString(err), fname, line);
    exit(1);
  }						  
}
