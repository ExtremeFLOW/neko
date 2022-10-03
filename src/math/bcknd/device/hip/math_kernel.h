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

/**
 * Device kernel for cmult
 */
template< typename T >
__global__ void cmult_kernel(T * __restrict__ a,
                             const T c,
                             const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = c * a[i];
  }
}

/**
 * Device kernel for cmult2
 */
template< typename T >
__global__ void cmult2_kernel(T * __restrict__ a,
                 T * __restrict__ b, 
                             const T c,
                             const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = c * b[i];
  }
}

/**
 * Device kernel for cadd
 */
template< typename T >
__global__ void cadd_kernel(T * __restrict__ a,
                            const T c,
                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c;
  }
}

/**
 * Device kernel for cmult
 */
template< typename T >
__global__ void cfill_kernel(T * __restrict__ a,
                             const T c,
                             const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = c;
  }
}

/**
 * Device kernel for add2
 */
template< typename T >
__global__ void add2_kernel(T * __restrict__ a,
                            const T * __restrict__ b,
                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + b[i];
  }
}

/**
 * Device kernel for add2s1
 */
template< typename T >
__global__ void add2s1_kernel(T * __restrict__ a,
                              const T * __restrict__ b,
                              const T c1,
                              const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = c1 * a[i] + b[i];
  }
}

/**
 * Device kernel for add2s2
 */
template< typename T >
__global__ void add2s2_many_kernel(T  * __restrict__  x,
                             const T ** p,
                             const T * alpha,
                 const int p_cur,
                             const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;


  for (int i = idx; i < n; i+= str) {
    T tmp = 0.0;
    for (int j = 0; j < p_cur; j ++) {
      tmp += p[j][i]*alpha[j];
    }
    x[i] += tmp;
  }
}

/**
 * Device kernel for add2s2
 */
template< typename T >
__global__ void add2s2_kernel(T * __restrict__ a,
                              const T * __restrict__ b,
                              const T c1,
                              const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c1 * b[i];
  }
}

/**
 * Device kernel for addsqr2s2
 */
template< typename T >
__global__ void addsqr2s2_kernel(T * __restrict__ a,
                                 const T * __restrict__ b,
                                 const T c1,
                                 const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c1 * (b[i] * b[i]);
  }
}

/**
 * Device kernel for add3s2
 */
template< typename T >
__global__ void add3s2_kernel(T * __restrict__ a,
                              const T * __restrict__ b,
                              const T * __restrict__ c,
                              const T c1,
                              const T c2,
                              const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = c1 * b[i] + c2 * c[i];
  }
}

/**
 * Device kernel for invcol1
 */
template< typename T >
__global__ void invcol1_kernel(T * __restrict__ a,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const T one = 1.0;

  for (int i = idx; i < n; i += str) {
    a[i] = one / a[i];
  }
}

/**
 * Device kernel for invcol2
 */
template< typename T >
__global__ void invcol2_kernel(T * __restrict__ a,
                               const T * __restrict__ b,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    a[i] = a[i] / b[i];
  }  
}

/** 
 * Device kernel for col2
 */
template< typename T >
__global__ void col2_kernel(T * __restrict__ a,
                            const T * __restrict__ b,
                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] * b[i];
  }  
}

/** 
 * Device kernel for col3
 */
template< typename T >
__global__ void col3_kernel(T * __restrict__ a,
                            const T * __restrict__ b,
                            const T * __restrict__ c,
                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] * c[i];
  }  
}

/** 
 * Device kernel for subcol3
 */
template< typename T >
__global__ void subcol3_kernel(T * __restrict__ a,
                               const T * __restrict__ b,
                               const T * __restrict__ c,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] - b[i] * c[i];
  }  
}

/** 
 * Device kernel for sub2
 */
template< typename T >
__global__ void sub2_kernel(T * __restrict__ a,
                            const T * __restrict__ b,
                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] - b[i];
  }  
}

/** 
 * Device kernel for sub3
 */
template< typename T >
__global__ void sub3_kernel(T * __restrict__ a,
                            const T * __restrict__ b,
                            const T * __restrict__ c,
                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] - c[i];
  }  
}

/**
 * Device kernel for addcol3
 */
template< typename T >
__global__ void addcol3_kernel(T * __restrict__ a,
                               const T * __restrict__ b,
                               const T * __restrict__ c,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + b[i] * c[i];
  }  
 
}

/**
 * Device kernel for addcol4
 */
template< typename T >
__global__ void addcol4_kernel(T * __restrict__ a,
                               const T * __restrict__ b,
                               const T * __restrict__ c,
                               const T * __restrict__ d,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + b[i] * c[i] * d[i];
  }  
  
}

/**
 * Warp shuffle reduction
 */
template< typename T>
__inline__ __device__ T reduce_warp(T val) {
  val += __shfl_down(val, 32);
  val += __shfl_down(val, 16);
  val += __shfl_down(val, 8);
  val += __shfl_down(val, 4);
  val += __shfl_down(val, 2);
  val += __shfl_down(val, 1);
  return val;
}

/**
 * Vector reduction kernel
 */
template< typename T >
__global__ void reduce_kernel(T * bufred, const int n) {
                
  T sum = 0;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  for (int i = idx; i<n ; i += str) 
  {
    sum += bufred[i];
  }

  __shared__ T shared[64];
  unsigned int lane = threadIdx.x % warpSize;
  unsigned int wid = threadIdx.x / warpSize;

  sum = reduce_warp<T>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    bufred[blockIdx.x] = sum;
}

/**
 * Reduction kernel for glsc3
 */

template< typename T >
__global__ void glsc3_reduce_kernel( T * bufred,
                                    const int n,
                                    const int j
                                   ) {
   __shared__ T buf[1024] ;
   const int idx = threadIdx.x;
   const int y= blockIdx.x;
   const int step = blockDim.x;

   buf[idx]=0;
   for (int i=idx ; i<n ; i+=step)
   {
     buf[idx] += bufred[i*j + y];
   }
   __syncthreads();

   int i = 512;
   while (i != 0)
   {
     if(threadIdx.x < i && (threadIdx.x + i) < n )
     {
        buf[threadIdx.x] += buf[threadIdx.x + i] ;
     }
     i = i>>1;
     __syncthreads();
   }

   bufred[y] = buf[0];
}

/**
 * Device kernel for glsc3
 */
template< typename T >
__global__ void glsc3_kernel(const T * a,
                             const T * b,
                             const T * c,
                             T * buf_h,
                             const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;
  
  __shared__ T shared[64];  
  T sum = 0.0;
  for (int i = idx; i < n; i+= str) {
    sum += a[i] * b[i] * c[i];
  }

  sum = reduce_warp<T>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    buf_h[blockIdx.x] = sum; 
}

/**
 * Device kernel for glsc3 many
 */
template< typename T >
__global__ void glsc3_many_kernel(const T * a,
                                  const T ** b,
                                  const T * c,
                                  T * buf_h,
                                  const int j,
                                  const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const int y = threadIdx.y;
  
  __shared__ T buf[1024];
  T tmp = 0;
  
  if(y < j){
    for (int i = idx; i < n; i+= str) {
      tmp += a[i] * b[threadIdx.y][i] * c[i];
    }
  }

  buf[threadIdx.x*blockDim.y+y] = tmp;
  __syncthreads();

  int i = blockDim.x>>1;
  while (i != 0) {
    if (threadIdx.x < i) {
      buf[threadIdx.x*blockDim.y +y] += buf[(threadIdx.x + i)*blockDim.y+y];
    }
    __syncthreads();
    i = i>>1;
  }
  
  if (threadIdx.x == 0) {
    if( y < j) {
      buf_h[j*blockIdx.x+y] = buf[y];
    }
  }
}

/**
 * Device kernel for glsc2
 */
template< typename T >
__global__ void glsc2_kernel(const T * a,
                             const T * b,
                             T * buf_h,
                             const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;
  
  __shared__ T shared[64];  
  T sum = 0.0;
  for (int i = idx; i < n; i+= str) {
    sum += a[i] * b[i];
  }

  sum = reduce_warp<T>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    buf_h[blockIdx.x] = sum;
}

/**
 * Device kernel for glsum
 */
template< typename T >
__global__ void glsum_kernel(const T * a,
                             T * buf_h,
                             const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;
  
  __shared__ T shared[64];
  T sum = 0;    
  for (int i = idx; i<n ; i += str) 
  {
    sum += a[i];
  }

  sum = reduce_warp<T>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    buf_h[blockIdx.x] = sum;

}
