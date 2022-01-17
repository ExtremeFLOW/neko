#ifndef __CUDA_CHECK_H
#define __CUDA_CHECK_H

void cuda_check(const char *fname, const int line, const cudaError_t err);

#define CUDA_CHECK(err) cuda_check(__FILE__, __LINE__, err)

#endif
