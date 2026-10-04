/*
 Copyright (c) 2026, The Neko Authors
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

/**
 * Device property queries used for hardware-aware work decomposition
 * (e.g. the element chunking in the dealiased advection), see
 * device_mp_count() and device_total_mem() in device.F90.
 */

#include <stdio.h>
#include <device/cuda/props.h>

extern "C" {

  /**
   * Number of streaming multiprocessors of the current device,
   * or 0 if the query fails.
   */
  int cuda_device_mp_count(void)
  {
    int dev = 0, count = 0;
    if (cudaGetDevice(&dev) != cudaSuccess) {
      (void) cudaGetLastError();
      return 0;
    }
    if (cudaDeviceGetAttribute(&count, cudaDevAttrMultiProcessorCount,
                               dev) != cudaSuccess) {
      (void) cudaGetLastError();
      return 0;
    }
    return count;
  }

  /**
   * Total memory of the current device in bytes, or 0 if the query fails.
   */
  size_t cuda_device_total_mem(void)
  {
    size_t mem_free = 0, mem_total = 0;
    if (cudaMemGetInfo(&mem_free, &mem_total) != cudaSuccess) {
      (void) cudaGetLastError();
      return 0;
    }
    return mem_total;
  }

}
