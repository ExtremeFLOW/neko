/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of the authors nor the names of its contributors may be
     used to endorse or promote products derived from this software without
     specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
*/

#include <device/hip/check.h>

using idp_stream_t = hipStream_t;

static inline void idp_zero(void *ptr, size_t bytes, idp_stream_t stream) {
  HIP_CHECK(hipMemsetAsync(ptr, 0, bytes, stream));
}

static inline void idp_check_launch() {
  HIP_CHECK(hipGetLastError());
}

#include "../euler_idp_kernel.h"
