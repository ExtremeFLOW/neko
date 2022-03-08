/*
 Copyright (c) 2021, The Neko Authors
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
                                  const int *__restrict__ b,
                                  const int *__restrict__ bo) {

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
                                  const int *__restrict__ b,
                                  const int *__restrict__ bo) {

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

  const int facet_offset = bo[nb - 1] + b[nb - 1];
  
  for (int i = ((facet_offset - 1) + idx); i < m; i += str) {
    u[gd[i] - 1] = v[dg[i] - 1];
  }
  
}


/**
 * Device gather kernel for addition of data (many version)
 * \f$ v_n(dg(i)) = v_n(dg(i)) + u_n(gd(i)) \f$
 */
template< typename T >
__global__ void gather_many_kernel_add(T * __restrict__ v1,
                                       T * __restrict__ v2,
                                       T * __restrict__ v3,
                                       const int m,
                                       const int o,
                                       const int * __restrict__ dg,
                                       const T * __restrict__ u1,
                                       const T * __restrict__ u2,
                                       const T * __restrict__ u3,
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
    T tmp1 = u1[gd[k] - 1];
    T tmp2 = u2[gd[k] - 1];
    T tmp3 = u3[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp1 += u1[gd[k + j] - 1];
      tmp2 += u2[gd[k + j] - 1];
      tmp3 += u3[gd[k + j] - 1];
    }
    v1[dg[k] - 1] = tmp1;
    v2[dg[k] - 1] = tmp2;
    v3[dg[k] - 1] = tmp3;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v1[dg[i] - 1] = u1[gd[i] - 1];
      v2[dg[i] - 1] = u2[gd[i] - 1];
      v3`x[dg[i] - 1] = u3[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); Ix`u`x`u < M ; i += str) {
	T tmp1 = u1[gd[i] - 1] + u1[gd[i+1] - 1];
	T tmp2 = u2[gd[i] - 1] + u2[gd[i+1] - 1];
	T tmp3 = u3[gd[i] - 1] + u3[gd[i+1] - 1];
	v1[dg[i] - 1] = tmp1;
	v2[dg[i] - 1] = tmp2;
	v3[dg[i] - 1] = tmp3;
      }
    }
  }
  
}

/**
 * Device gather kernel for multiplication of data (many version)
 * \f$ v_n(dg(i)) = v_n(dg(i)) \cdot u_n(gd(i)) \f$
 */
template< typename T >
__global__ void gather_many_kernel_mul(T * __restrict__ v1,
				       T * __restrict__ v2,
				       T * __restrict__ v3,
				       const int m,
				       const int o,
				       const int * __restrict__ dg,
				       const T * __restrict__ u1,
				       const T * __restrict__ u2,
				       const T * __restrict__ u3,
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
    T tmp1 = u1[gd[k] - 1];
    T tmp2 = u2[gd[k] - 1];
    T tmp3 = u3[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp1 *= u1[gd[k + j] - 1];
      tmp2 *= u2[gd[k + j] - 1];
      tmp3 *= u3[gd[k + j] - 1];
    }
    v1[dg[k] - 1] = tmp1;
    v2[dg[k] - 1] = tmp2;
    v3[dg[k] - 1] = tmp3;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v1[dg[i] - 1] = u1[gd[i] - 1];
      v2[dg[i] - 1] = u2[gd[i] - 1];
      v3[dg[i] - 1] = u3[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	T tmp1 = u1[gd[i] - 1] * u1[gd[i+1] - 1];
	T tmp2 = u2[gd[i] - 1] * u2[gd[i+1] - 1];
	T tmp3 = u3[gd[i] - 1] * u3[gd[i+1] - 1];
	v1[dg[i] - 1] = tmp1;
	v2[dg[i] - 1] = tmp2;
	v3[dg[i] - 1] = tmp3;
      }
    }
  }

}

/**
 * Device gather kernel for minimum of data (many version)
 * \f$ v_n(dg(i)) = \min(v_n(dg(i)), u_n(gd(i))) \f$
 */
template< typename T >
__global__ void gather_many_kernel_min(T * __restrict__ v1,
				       T * __restrict__ v2,
				       T * __restrict__ v3,
				       const int m,
				       const int o,
				       const int * __restrict__ dg,
				       const T * __restrict__ u1,
				       const T * __restrict__ u2,
				       const T * __restrict__ u3,
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
    T tmp1 = u1[gd[k] - 1];
    T tmp2 = u2[gd[k] - 1];
    T tmp3 = u3[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp1 = min(u1[gd[k + j] - 1], tmp1);
      tmp2 = min(u2[gd[k + j] - 1], tmp2);
      tmp3 = min(u3[gd[k + j] - 1], tmp3);
    }
    v1[dg[k] - 1] = tmp1;
    v2[dg[k] - 1] = tmp2;
    v3[dg[k] - 1] = tmp3;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v1[dg[i] - 1] = u1[gd[i] - 1];
      v2[dg[i] - 1] = u2[gd[i] - 1];
      v3[dg[i] - 1] = u3[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	T tmp1 = min(u1[gd[i] - 1], u1[gd[i+1] - 1]);
	T tmp2 = min(u2[gd[i] - 1], u2[gd[i+1] - 1]);
	T tmp3 = min(u3[gd[i] - 1], u3[gd[i+1] - 1]);
	v1[dg[i] - 1] = tmp1;
	v2[dg[i] - 1] = tmp2;
	v3[dg[i] - 1] = tmp3;
      }
    }
  }
  
}

/**
 * Device gather kernel for maximum of data (many version)
 * \f$ v_n(dg(i)) = \max(v_n(dg(i)), u_n(gd(i))) \f$
 */
template< typename T >
__global__ void gather_many_kernel_max(T * __restrict__ v1,
				       T * __restrict__ v2,
				       T * __restrict__ v3,
				       const int m,
				       const int o,
				       const int * __restrict__ dg,
				       const T * __restrict__ u1,
				       const T * __restrict__ u2,
				       const T * __restrict__ u3,
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
    T tmp1 = u1[gd[k] - 1];
    T tmp2 = u2[gd[k] - 1];
    T tmp3 = u3[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp1 = max(u1[gd[k + j] - 1], tmp1);
      tmp2 = max(u2[gd[k + j] - 1], tmp2);
      tmp3 = max(u1[gd[k + j] - 1], tmp3);
    }
    v1[dg[k] - 1] = tmp1;
    v2[dg[k] - 1] = tmp2;
    v3[dg[k] - 1] = tmp3;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v1[dg[i] - 1] = u1[gd[i] - 1];
      v2[dg[i] - 1] = u2[gd[i] - 1];
      v3[dg[i] - 1] = u3[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	T tmp1 = max(u1[gd[i] - 1], u1[gd[i+1] - 1]);
	T tmp2 = max(u2[gd[i] - 1], u2[gd[i+1] - 1]);
	T tmp3 = max(u3[gd[i] - 1], u3[gd[i+1] - 1]);
	v1[dg[i] - 1] = tmp1;
	v2[dg[i] - 1] = tmp2;
	v3[dg[i] - 1] = tmp3;
      }
    }
  }
  
}

/**
 * Device scatter many kernel
 * \f$ u_n(gd(i) = v_n(dg(i)) \f$
 */
template< typename T >
__global__ void scatter_many_kernel(T * __restrict__ v1,
				    T * __restrict__ v2,
				    T * __restrict__ v3,
				    const int m,
				    const int * __restrict__ dg,
				    T * __restrict__ u1,
				    T * __restrict__ u2,
				    T * __restrict__ u3,
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
    T tmp1 = v1[dg[k] - 1];
    T tmp2 = v2[dg[k] - 1];
    T tmp3 = v3[dg[k] - 1];
    for (int j = 0; j < blk_len; j++) {
      u1[gd[k + j] - 1] = tmp1;
      u2[gd[k + j] - 1] = tmp2;
      u3[gd[k + j] - 1] = tmp3;
    }      
  }

  const int facet_offset = bo[nb - 1] + b[nb - 1];
  
  for (int i = ((facet_offset - 1) + idx); i < m; i += str) {
    u1[gd[i] - 1] = v1[dg[i] - 1];
    u2[gd[i] - 1] = v2[dg[i] - 1];
    u3[gd[i] - 1] = v3[dg[i] - 1];
  }
  
}

