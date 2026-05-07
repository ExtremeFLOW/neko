/*
 Copyright (c) 2021-2022, The Neko Authors
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

#ifndef __GS_GS_KERNELS__
#define __GS_GS_KERNELS__

/**
 * Device gather kernel for addition of data
 * \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
 */
template< typename T >
__global__ void gather_kernel_add(T * __restrict__ v,
                                  const int m,
                                  const int o,
                                  const int * __restrict__ dg,
                                  const T * __restrict__ u,
                                  const int n,
                                  const int * __restrict__ gd,
                                  const int nb,
                                  const int * __restrict__ b,
                                  const int * __restrict__ bo) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    T tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp += u[gd[k + j] - 1];
    }
    v[dg[k] - 1] = tmp;
  }

  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
        T tmp = u[gd[i] - 1] + u[gd[i+1] - 1];
        v[dg[i] - 1] = tmp;
      }
    }
  }

}

/**
 * Device gather kernel for multiplication of data
 * \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
 */
template< typename T >
__global__ void gather_kernel_mul(T * __restrict__ v,
                                  const int m,
                                  const int o,
                                  const int * __restrict__ dg,
                                  const T * __restrict__ u,
                                  const int n,
                                  const int * __restrict__ gd,
                                  const int nb,
                                  const int * __restrict__ b,
                                  const int * __restrict__ bo) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    T tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp *= u[gd[k + j] - 1];
    }
    v[dg[k] - 1] = tmp;
  }

  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
        T tmp = u[gd[i] - 1] * u[gd[i+1] - 1];
        v[dg[i] - 1] = tmp;
      }
    }
  }

}

/**
 * Device gather kernel for minimum of data
 * \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
 */
template< typename T >
__global__ void gather_kernel_min(T * __restrict__ v,
                                  const int m,
                                  const int o,
                                  const int * __restrict__ dg,
                                  const T * __restrict__ u,
                                  const int n,
                                  const int * __restrict__ gd,
                                  const int nb,
                                  const int * __restrict__ b,
                                  const int * __restrict__ bo) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    T tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp = min(u[gd[k + j] - 1], tmp);
    }
    v[dg[k] - 1] = tmp;
  }

  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
        T tmp = min(u[gd[i] - 1], u[gd[i+1] - 1]);
        v[dg[i] - 1] = tmp;
      }
    }
  }

}

/**
 * Device gather kernel for maximum of data
 * \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
 */
template< typename T >
__global__ void gather_kernel_max(T * __restrict__ v,
                                  const int m,
                                  const int o,
                                  const int * __restrict__ dg,
                                  const T * __restrict__ u,
                                  const int n,
                                  const int * __restrict__ gd,
                                  const int nb,
                                  const int * __restrict__ b,
                                  const int * __restrict__ bo) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    T tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp = max(u[gd[k + j] - 1], tmp);
    }
    v[dg[k] - 1] = tmp;
  }

  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
        T tmp = max(u[gd[i] - 1], u[gd[i+1] - 1]);
        v[dg[i] - 1] = tmp;
      }
    }
  }

}

/**
 * Device scatter kernel
 * \f$ u(gd(i) = v(dg(i)) \f$
 */
template< typename T >
__global__ void scatter_kernel(T * __restrict__ v,
                               const int m,
                               const int * __restrict__ dg,
                               T * __restrict__ u,
                               const int n,
                               const int * __restrict__ gd,
                               const int nb,
                               const int *__restrict__ b,
                               const int *__restrict__ bo) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    T tmp = v[dg[k] - 1];
    for (int j = 0; j < blk_len; j++) {
      u[gd[k + j] - 1] = tmp;
    }
  }

  int facet_offset = 0;
  if ( nb > 0) {
    facet_offset = bo[nb - 1] + b[nb - 1];
  }
  

  for (int i = (facet_offset + idx); i < m; i += str) {
    u[gd[i] - 1] = v[dg[i] - 1];
  }

}

template< typename T >
__global__ void gs_pack_kernel(const T * __restrict__ u,
                               T * __restrict__ buf,
                               const int32_t * __restrict__ dof,
                               const int n) {

  const int j = threadIdx.x + blockDim.x * blockIdx.x;

  if (j >= n)
    return;

  buf[j] = u[dof[j]-1];
}


template< typename T >
__global__ void gs_unpack_add_kernel(T * __restrict__ u,
                                     const T * __restrict__ buf,
                                     const int32_t * __restrict__ dof,
                                     const int n) {

  const int j = threadIdx.x + blockDim.x * blockIdx.x;

  if (j >= n)
    return;

  const int32_t idx = dof[j];
  const T val = buf[j];
  if (idx < 0) {
    atomicAdd(&u[-idx-1], val);
  } else {
    u[idx-1] += val;
  }
}

template<typename T>
__device__ T atomicMinFloat(T* address, T val);

template<>
__device__ float atomicMinFloat<float>(float* address, float val) {
    float old;
    old = !signbit(val) ? __int_as_float(atomicMin((int*)address, __float_as_int(val))) :
        __uint_as_float(atomicMax((unsigned int*)address, __float_as_uint(val)));

    return old;
}

template<>
__device__ double atomicMinFloat<double>(double* address, double val) {
    double old;
    old = !signbit(val) ? __longlong_as_double(atomicMin((unsigned long long*)address, 
                                            __double_as_longlong(val))) :
          __longlong_as_double(atomicMax((unsigned long long*)address, 
                                            __double_as_longlong(val)));
    return old;
}

template< typename T >
__global__ void gs_unpack_min_kernel(T * __restrict__ u,
                                     const T * __restrict__ buf,
                                     const int32_t * __restrict__ dof,
                                     const int n) {

  const int j = threadIdx.x + blockDim.x * blockIdx.x;

  if (j >= n)
    return;

  const int32_t idx = dof[j];
  const T val = buf[j];

  if (idx < 0) {
    // Use atomicMin for shared nodal points
    atomicMinFloat(&u[-idx-1], val);
  } else {
    // Directly compute min for nodal points on edges
    u[idx-1] = min(u[idx-1], val);
  }
}

template<typename T>
__device__ T atomicMaxFloat(T* address, T val);

template<>
__device__ float atomicMaxFloat<float>(float* address, float val) {
    float old;
    old = !signbit(val) ? __int_as_float(atomicMax((int*)address, __float_as_int(val))) :
        __uint_as_float(atomicMin((unsigned int*)address, __float_as_uint(val)));
    return old;
}

template<>
__device__ double atomicMaxFloat<double>(double* address, double val) {
    double old;
    old = !signbit(val) ? __longlong_as_double(atomicMax((unsigned long long*)address, 
                                                __double_as_longlong(val))) :
          __longlong_as_double(atomicMin((unsigned long long*)address, 
                                                __double_as_longlong(val)));
    return old;
}

template< typename T >
__global__ void gs_unpack_max_kernel(T * __restrict__ u,
                                     const T * __restrict__ buf,
                                     const int32_t * __restrict__ dof,
                                     const int n) {

  const int j = threadIdx.x + blockDim.x * blockIdx.x;

  if (j >= n)
    return;

  const int32_t idx = dof[j];
  const T val = buf[j];

  if (idx < 0) {
    // Use atomicMax for shared nodal points
    atomicMaxFloat(&u[-idx-1], val);
  } else {
    // Directly compute min for nodal points on edges
    u[idx-1] = max(u[idx-1], val);
  }
}

#endif // __GS_GS_KERNELS__