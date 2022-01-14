#ifndef __HIP_CHECK_H
#define __HIP_CHECK_H

#include <hip/hip_runtime.h>

void hip_check(const char *fname, const int line, const hipError_t err);

#define HIP_CHECK(err) hip_check(__FILE__, __LINE__, err)

#endif
