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

#include <cmath>
#include <cfloat>
#include <device/device_config.h>

namespace {

template<typename T> __device__ __forceinline__ T real_epsilon();
template<typename T> __device__ __forceinline__ T real_tiny();
template<typename T> __device__ __forceinline__ T real_huge();

template<> __device__ __forceinline__ float real_epsilon<float>() {
  return FLT_EPSILON;
}

template<> __device__ __forceinline__ double real_epsilon<double>() {
  return DBL_EPSILON;
}

template<> __device__ __forceinline__ float real_tiny<float>() {
  return FLT_MIN;
}

template<> __device__ __forceinline__ double real_tiny<double>() {
  return DBL_MIN;
}

template<> __device__ __forceinline__ float real_huge<float>() {
  return FLT_MAX;
}

template<> __device__ __forceinline__ double real_huge<double>() {
  return DBL_MAX;
}

template<typename T> __device__ __forceinline__ T atomic_min_real(T *p, T v);
template<typename T> __device__ __forceinline__ T atomic_max_real(T *p, T v);

template<> __device__ __forceinline__ float atomic_min_real(float *p, float v) {
  int *address = reinterpret_cast<int *>(p);
  int old = *address;
  while (v < __int_as_float(old)) {
    const int assumed = old;
    old = atomicCAS(address, assumed, __float_as_int(v));
    if (old == assumed) break;
  }
  return __int_as_float(old);
}

template<> __device__ __forceinline__ float atomic_max_real(float *p, float v) {
  int *address = reinterpret_cast<int *>(p);
  int old = *address;
  while (v > __int_as_float(old)) {
    const int assumed = old;
    old = atomicCAS(address, assumed, __float_as_int(v));
    if (old == assumed) break;
  }
  return __int_as_float(old);
}

template<> __device__ __forceinline__ double atomic_min_real(double *p,
                                                               double v) {
  unsigned long long *address = reinterpret_cast<unsigned long long *>(p);
  unsigned long long old = *address;
  while (v < __longlong_as_double(old)) {
    const unsigned long long assumed = old;
    old = atomicCAS(address, assumed, __double_as_longlong(v));
    if (old == assumed) break;
  }
  return __longlong_as_double(old);
}

template<> __device__ __forceinline__ double atomic_max_real(double *p,
                                                               double v) {
  unsigned long long *address = reinterpret_cast<unsigned long long *>(p);
  unsigned long long old = *address;
  while (v > __longlong_as_double(old)) {
    const unsigned long long assumed = old;
    old = atomicCAS(address, assumed, __double_as_longlong(v));
    if (old == assumed) break;
  }
  return __longlong_as_double(old);
}


#define IDP_GLOBAL
#define IDP_PRIVATE
#define IDP_SHARED __shared__
#define IDP_VALUE(type,name) type name
#define IDP_GRID
#define IDP_GROUP blockIdx.x
#define IDP_WIDTH blockDim.x
#define IDP_GROUPS gridDim.x
#define IDP_LANE threadIdx.x
#define IDP_BARRIER() __syncthreads()
#define IDP_KERNEL __global__
#define IDP_INLINE __device__ inline
#define IDP_EPS real_epsilon<real>()
#define IDP_TINY real_tiny<real>()
#define IDP_HUGE real_huge<real>()
#define idp_add atomicAdd
#define idp_minimum atomic_min_real
#define idp_maximum atomic_max_real
#include "euler_idp_kernels.cl"

} // namespace

#define RBUF(p) ((real *)(p))
#define IBUF(p) ((int *)(p))
#define RVAL(x) (x)
#define IVAL(x) (x)
#define IDP_ZERO(p,n) idp_zero(p, sizeof(real)*(n), (idp_stream_t)glb_cmd_queue)
#define IDP_LAUNCH(kernel,groups,width,...) do { \
  kernel<<<groups,width,0,(idp_stream_t)glb_cmd_queue>>>(__VA_ARGS__); \
  idp_check_launch(); \
} while (0)
extern "C" {
#include "euler_idp_launch.h"
}
