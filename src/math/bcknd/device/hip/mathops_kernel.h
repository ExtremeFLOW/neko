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
 * Device kernel for opchsign
 */
template< typename T>
__global__ void opchsign_kernel(T * __restrict__ a1, 
                                T * __restrict__ a2, 
                                T * __restrict__ a3,
                                const int gdim,
                                const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = -a1[i];
      a2[i] = -a2[i];
      a3[i] = -a3[i];
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = -a1[i];
      a2[i] = -a2[i];
    }
  }

}

/**
 * Device kernel for opcolv
 */
template< typename T>
__global__ void opcolv_kernel(T * __restrict__ a1, 
                              T * __restrict__ a2, 
                              T * __restrict__ a3,
                              const T * __restrict__ c,
                              const int gdim,
                              const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] * c[i];
      a2[i] = a2[i] * c[i];
      a3[i] = a3[i] * c[i];
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] * c[i];
      a2[i] = a2[i] * c[i];
    }
  }

}

/**
 * Device kernel for opcolv3c
 */
template< typename T>
__global__ void opcolv3c_kernel(T * __restrict__ a1, 
                                T * __restrict__ a2, 
                                T * __restrict__ a3,
                                const T * __restrict__ b1, 
                                const T * __restrict__ b2, 
                                const T * __restrict__ b3,
                                const T * __restrict__ c,
                                const T d,
                                const int gdim,
                                const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = b1[i] * c[i] * d;
      a2[i] = b2[i] * c[i] * d;
      a3[i] = b3[i] * c[i] * d;
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = b1[i] * c[i] * d;
      a2[i] = b2[i] * c[i] * d;
    }
  }

}

/**
 * Device kernel for opadd2cm
 */
template< typename T>
__global__ void opadd2cm_kernel(T * __restrict__ a1, 
                                T * __restrict__ a2, 
                                T * __restrict__ a3,
                                const T * __restrict__ b1, 
                                const T * __restrict__ b2, 
                                const T * __restrict__ b3,
                                const T c,
                                const int gdim,
                                const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] + b1[i] * c;
      a2[i] = a2[i] + b2[i] * c;
      a3[i] = a3[i] + b3[i] * c;
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] + b1[i] * c;
      a2[i] = a2[i] + b2[i] * c;
    }
  }

}


/**
 * Device kernel for opadd2col
 */
template< typename T>
__global__ void opadd2col_kernel(T * __restrict__ a1, 
                                 T * __restrict__ a2, 
                                 T * __restrict__ a3,
                                 const T * __restrict__ b1, 
                                 const T * __restrict__ b2, 
                                 const T * __restrict__ b3,
                                 const T * __restrict__ c,
                                 const int gdim,
                                 const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] + b1[i] * c[i];
      a2[i] = a2[i] + b2[i] * c[i];
      a3[i] = a3[i] + b3[i] * c[i];
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] + b1[i] * c[i];
      a2[i] = a2[i] + b2[i] * c[i];
    }
  }

}

