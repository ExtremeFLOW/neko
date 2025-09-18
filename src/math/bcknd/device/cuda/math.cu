/*
 Copyright (c) 2021-2025, The Neko Authors
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

#include "math_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_NVSHMEM
#include <nvshmem.h>
#include <nvshmemx.h>
#endif

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

#ifdef HAVE_NCCL
#include <math/bcknd/device/device_nccl_reduce.h>
#include <math/bcknd/device/device_nccl_op.h>
#endif

  /** Fortran wrapper for copy
   * Copy a vector \f$ a = b \f$
   */
  void cuda_copy(void *a, void *b, int *n, cudaStream_t strm) {
    CUDA_CHECK(cudaMemcpyAsync(a, b, (*n) * sizeof(real),
                               cudaMemcpyDeviceToDevice, strm));
  }

  /** Fortran wrapper for masked copy
   * Copy a vector \f$ a(mask) = b(mask) \f$
   */
  void cuda_masked_copy(void *a, void *b, void *mask,
                        int *n, int *m, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    masked_copy_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real*) b,(int*) mask, *n, *m);
    CUDA_CHECK(cudaGetLastError());

  }

  /** Fortran wrapper for masked gather copy
   * Copy a vector \f$ a(i) = b(mask(i)) \f$
   */
  void cuda_masked_gather_copy(void *a, void *b, void *mask,
                               int *n, int *m, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    masked_gather_copy_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real*) b,(int*) mask, *n, *m);
    CUDA_CHECK(cudaGetLastError());

  }

  /** Fortran wrapper for masked atomic reduction
   * update a vector \f$ a += b(mask) \f$ where mask is not unique
   */
  void cuda_masked_atomic_reduction(void *a, void *b, void *mask,
                                    int *n, int *m, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    masked_atomic_reduction_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (int *) mask, *n, *m);
    CUDA_CHECK(cudaGetLastError());

  }
  /** Fortran wrapper for masked scatter copy
   * Copy a vector \f$ a(mask(i)) = b(i) \f$
   */
  void cuda_masked_scatter_copy(void *a, void *b, void *mask,
                                int *n, int *m, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    masked_scatter_copy_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real*) b,(int*) mask, *n, *m);
    CUDA_CHECK(cudaGetLastError());
  }


  /** Fortran wrapper for cfill_mask
   * Fill a scalar to vector \f$ a_i = s, for i \in mask \f$
   */
  void cuda_cfill_mask(void* a, real* c, int* size, int* mask, int* mask_size,
                       cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*mask_size) + 1024 - 1) / 1024, 1, 1);

    cfill_mask_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real*)a, *c, *size, mask, *mask_size);
    CUDA_CHECK(cudaGetLastError());
  }

  /** Fortran wrapper for rzero
   * Zero a real vector
   */
  void cuda_rzero(void *a, int *n, cudaStream_t strm) {
    CUDA_CHECK(cudaMemsetAsync(a, 0, (*n) * sizeof(real), strm));
  }

  /** Fortran wrapper for cmult
   * Multiplication by constant c \f$ a = c \cdot a \f$
   */
  void cuda_cmult(void *a, real *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    cmult_kernel<real><<<nblcks, nthrds, 0, strm>>>((real *) a, *c, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /** Fortran wrapper for cmult2
   * Multiplication by constant c \f$ a = c \cdot b \f$
   */
  void cuda_cmult2(void *a, void *b, real *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    cmult2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, *c, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /** Fortran wrapper for cdiv
   * Division of constant c by array \f$ a = c / a \f$
   */
  void cuda_cdiv(void *a, real *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    cdiv_kernel<real><<<nblcks, nthrds, 0, strm>>>((real *) a, *c, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /** Fortran wrapper for cdiv2
   * Division of constant c by array \f$ a = c / b \f$
   */
  void cuda_cdiv2(void *a, void *b, real *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    cdiv2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, *c, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /** Fortran wrapper for cadd
   * Add a scalar to vector \f$ a_i = a_i + c \f$
   */
  void cuda_radd(void *a, real *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    cadd_kernel<real><<<nblcks, nthrds, 0, strm>>>((real *) a, *c, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for cadd2
   * Add a scalar to vector \f$ a_i = b_i + c \f$
   */
  void cuda_cadd2(void *a, void *b, real *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    cadd2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, *c, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /** Fortran wrapper for cfill
   * Set all elements to a constant c \f$ a = c \f$
   */
  void cuda_cfill(void *a, real *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    if (*n > 0){
      cfill_kernel<real><<<nblcks, nthrds, 0, strm>>>((real *) a, *c, *n);
      CUDA_CHECK(cudaGetLastError());
    }

  }

  /**
   * Fortran wrapper for add2
   * Vector addition \f$ a = a + b \f$
   */
  void cuda_add2(void *a, void *b, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for add3
   * Vector addition \f$ a = b + c \f$
   */
  void cuda_add3(void *a, void *b, void *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add3_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for add4
   * Vector addition \f$ a = b + c + d \f$
   */
  void cuda_add4(void *a, void *b, void *c, void *d, int *n,
                 cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add4_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, (real *) d, *n);
    CUDA_CHECK(cudaGetLastError());

  }
  /**
   * Fortran wrapper for add2s1
   * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
   * (multiplication on first argument)
   */
  void cuda_add2s1(void *a, void *b, real *c1, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add2s1_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, *c1, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for add2s2
   * Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
   * (multiplication on second argument)
   */
  void cuda_add2s2(void *a, void *b, real *c1, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add2s2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, *c1, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for add2s2
   * Vector addition with scalar multiplication
   * \f$ x = x + c_1 p1 + c_2p2 + ... + c_jpj \f$
   * (multiplication on second argument)
   */
  void cuda_add2s2_many(void *x, void **p, void *alpha, int *j, int *n,
                        cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add2s2_many_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) x, (const real **) p, (real *) alpha, *j, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for addsqr2s2
   * Vector addition with scalar multiplication \f$ a = a + c_1 (b * b) \f$
   * (multiplication on second argument)
   */
  void cuda_addsqr2s2(void *a, void *b, real *c1, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    addsqr2s2_kernel<real><<<nblcks, nthrds, 0,strm>>>
      ((real *) a, (real *) b, *c1, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for add3s2
   * Vector addition with scalar multiplication \f$ a = c_1 b + c_2 c \f$
   * (multiplication on second argument)
   */
  void cuda_add3s2(void *a, void *b, void *c, real *c1, real *c2, int *n,
                   cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add3s2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, *c1, *c2, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for add4s3
   * Vector addition with scalar multiplication \f$ a = c_1 b + c_2 c + c_3 d\f$
   * (multiplication on second argument)
   */
  void cuda_add4s3(void *a, void *b, void *c, void *d, real *c1, real *c2,
                   real *c3, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add4s3_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, (real *) d, *c1, *c2, *c3, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for add5s4
   * Vector addition with scalar multiplication \f$ a = a + c_1 b + c_2 c + c_3 d + c_4 e\f$
   * (multiplication on second argument)
   */
  void cuda_add5s4(void *a, void *b, void *c, void *d, void *e, real *c1,
                   real *c2, real *c3, real *c4, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    add5s4_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, (real *) d, (real *) e,
       *c1, *c2, *c3, *c4, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper for invcol1
   * Invert a vector \f$ a = 1 / a \f$
   */
  void cuda_invcol1(void *a, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    invcol1_kernel<real><<<nblcks, nthrds, 0, strm>>>((real *) a, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for invcol2
   * Vector division \f$ a = a / b \f$
   */
  void cuda_invcol2(void *a, void *b, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    invcol2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for invcol3
   * Vector division \f$ a = b / c \f$
   */
  void cuda_invcol3(void *a, void *b, void *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    invcol3_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b,  (real *) c, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for col2
   * Vector multiplication with 2 vectors \f$ a = a \cdot b \f$
   */
  void cuda_col2(void *a, void *b, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    col2_kernel<real><<<nblcks, nthrds, 0, strm>>>((real *) a, (real *) b, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for col3
   * Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
   */
  void cuda_col3(void *a, void *b, void *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    col3_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for subcol3
   * Vector multiplication with 3 vectors \f$ a = a - b \cdot c \f$
   */
  void cuda_subcol3(void *a, void *b, void *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    subcol3_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, *n);
    CUDA_CHECK(cudaGetLastError());
  }


  /**
   * Fortran wrapper for sub2
   * Vector subtraction \f$ a = a - b \f$
   */
  void cuda_sub2(void *a, void *b, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    sub2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for sub3
   * Vector subtraction \f$ a = b - c \f$
   */
  void cuda_sub3(void *a, void *b, void *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    sub3_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for addcol3
   * \f$ a = a + b * c \f$
   */
  void cuda_addcol3(void *a, void *b, void *c, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    addcol3_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for addcol4
   * \f$ a = a + b * c * d\f$
   */
  void cuda_addcol4(void *a, void *b, void *c, void *d, int *n,
                    cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    addcol4_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, (real *) d, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for addcol3s2
   * \f$ a = a + s(b * c) \f$
   */
  void cuda_addcol3s2(void *a, void *b, void *c, real *s, int *n,
                      cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    addcol3s2_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) a, (real *) b, (real *) c, *s, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for vdot3
   * \f$ dot = u \cdot v \f$
   */
  void cuda_vdot3(void *dot, void *u1, void *u2, void *u3,
                  void *v1, void *v2, void *v3, int *n,
                  cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    vdot3_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) dot, (real *) u1, (real *) u2, (real *) u3,
       (real *) v1, (real *) v2, (real *) v3, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fortran wrapper for vcross
   * \f$ u = v \times w \f$
   */
  void cuda_vcross(void *u1, void *u2, void *u3,
                  void *v1, void *v2, void *v3,
                  void *w1, void *w2, void *w3,
                   int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    vcross_kernel<real><<<nblcks, nthrds, 0, strm>>>
      ((real *) u1, (real *) u2, (real *) u3,
       (real *) v1, (real *) v2, (real *) v3,
       (real *) w1, (real *) w2, (real *) w3, *n);
    CUDA_CHECK(cudaGetLastError());
  }


  /*
   * Reduction buffer
   */
  int red_s = 0;
  real * bufred = NULL;
  void * bufred_d = NULL;

  /**
   *Checks and allocates a buffer of size nb*sizeof(real) for reductions
  */
  void cuda_redbuf_check_alloc(int nb) {
    if ( nb >= red_s) {
      red_s = nb+1;
      if (bufred != NULL) {
        CUDA_CHECK(cudaFreeHost(bufred));
#ifdef HAVE_NVSHMEM
        nvshmem_free(bufred_d);
#else
        CUDA_CHECK(cudaFree(bufred_d));
#endif
      }
      CUDA_CHECK(cudaMallocHost(&bufred,red_s*sizeof(real)));
#ifdef HAVE_NVSHMEM
      bufred_d = (real *) nvshmem_malloc(red_s*sizeof(real));
#else
      CUDA_CHECK(cudaMalloc(&bufred_d, red_s*sizeof(real)));
#endif
    }
  }

  /**
   * Global additive reduction
   */
  void cuda_global_reduce_add(real * bufred, void * bufred_d,
                              int n, const cudaStream_t stream) {

#ifdef HAVE_NCCL
    device_nccl_allreduce(bufred_d, bufred_d, n, sizeof(real),
                          DEVICE_NCCL_SUM, stream);
    CUDA_CHECK(cudaMemcpyAsync(bufred, bufred_d, sizeof(real)*n,
                               cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
#elif HAVE_NVSHMEM
    if (sizeof(real) == sizeof(float)) {
      nvshmemx_float_sum_reduce_on_stream(NVSHMEM_TEAM_WORLD,
                                           (float *) bufred_d,
                                           (float *) bufred_d, n, stream);
    }
    else if (sizeof(real) == sizeof(double)) {
      nvshmemx_double_sum_reduce_on_stream(NVSHMEM_TEAM_WORLD,
                                           (double *) bufred_d,
                                           (double *) bufred_d, n, stream);

    }
    CUDA_CHECK(cudaMemcpyAsync(bufred, bufred_d,
                               sizeof(real)*n, cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
#elif HAVE_DEVICE_MPI
    cudaStreamSynchronize(stream);
    device_mpi_allreduce(bufred_d, bufred, n, sizeof(real), DEVICE_MPI_SUM);
#else
    CUDA_CHECK(cudaMemcpyAsync(bufred, bufred_d, n*sizeof(real),
                               cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
#endif
  }



  /**
   * Fortran wrapper vlsc3
   * Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
   */
  real cuda_vlsc3(void *u, void *v, void *w, int *n, cudaStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    cuda_redbuf_check_alloc(nb);

    glsc3_kernel<real><<<nblcks, nthrds, 0, stream>>>
      ((real *) u, (real *) v, (real *) w, (real *) bufred_d, *n);
    CUDA_CHECK(cudaGetLastError());
    reduce_kernel<real><<<1, 1024, 0, stream>>> ((real *) bufred_d, nb);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemcpyAsync(bufred, bufred_d, sizeof(real),
                               cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);

    return bufred[0];
  }




  /**
   * Fortran wrapper glsc3
   * Weighted inner product \f$ a^T b c \f$
   */
  real cuda_glsc3(void *a, void *b, void *c, int *n, cudaStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    cuda_redbuf_check_alloc(nb);

    if ( *n > 0) {
      glsc3_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((real *) a, (real *) b, (real *) c, (real *) bufred_d, *n);
      CUDA_CHECK(cudaGetLastError());
      reduce_kernel<real><<<1, 1024, 0, stream>>> ((real *) bufred_d, nb);
      CUDA_CHECK(cudaGetLastError());
    }
    else {
      cuda_rzero(bufred_d, &red_s, stream);
    }
    cuda_global_reduce_add(bufred, bufred_d, 1, stream);

    return bufred[0];
  }

  /**
   * Fortran wrapper for doing an reduction to an array
   * Weighted inner product \f$ w^T v(n,1:j) c \f$
   */
  void cuda_glsc3_many(real *h, void * w, void *v,void *mult, int *j, int *n,
                       cudaStream_t stream){
    int pow2 = 1;
    while(pow2 < (*j)){
      pow2 = 2*pow2;
    }
    const int nt = 1024/pow2;
    const dim3 nthrds(pow2, nt, 1);
    const dim3 nblcks(((*n)+nt - 1)/nt, 1, 1);
    const int nb = ((*n) + nt - 1)/nt;

    cuda_redbuf_check_alloc((*j)*nb);

    if ( *n > 0) {
      glsc3_many_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((const real *) w, (const real **) v,
         (const real *)mult, (real *)bufred_d, *j, *n);
      CUDA_CHECK(cudaGetLastError());
      glsc3_reduce_kernel<real>
        <<<(*j), 1024, 0, stream>>>((real *) bufred_d, nb, *j);
      CUDA_CHECK(cudaGetLastError());
    }
    else {
      cuda_rzero(bufred_d, &red_s, stream);
    }
    cuda_global_reduce_add(h, bufred_d, (*j), stream);
  }

  /**
   * Fortran wrapper glsc2
   * Weighted inner product \f$ a^T b c \f$
   */
  real cuda_glsc2(void *a, void *b, int *n, cudaStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    cuda_redbuf_check_alloc(nb);

    if ( *n > 0) {
      glsc2_kernel<real>
        <<<nblcks, nthrds, 0, stream>>>((real *) a,
                                        (real *) b,
                                        (real *) bufred_d, *n);
      CUDA_CHECK(cudaGetLastError());
      reduce_kernel<real><<<1, 1024, 0, stream>>> ((real *) bufred_d, nb);
      CUDA_CHECK(cudaGetLastError());
    }
    else {
      cuda_rzero(bufred_d, &red_s, stream);
    }
    cuda_global_reduce_add(bufred, bufred_d, 1, stream);

    return bufred[0];
  }

  /**
   * Fortran wrapper glsubnorm
   * Squared Norm of difference \f$ \| a - b \|_2^2 \f$
   */
  real cuda_glsubnorm2(void *a, void *b, int *n, cudaStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    cuda_redbuf_check_alloc(nb);

    if ( *n > 0) {
      glsubnorm2_kernel<real>
        <<<nblcks, nthrds, 0, stream>>>((real *) a,
                                        (real *) b,
                                        (real *) bufred_d, *n);
      CUDA_CHECK(cudaGetLastError());
      reduce_kernel<real><<<1, 1024, 0, stream>>> ((real *) bufred_d, nb);
      CUDA_CHECK(cudaGetLastError());
    }
    else {
      cuda_rzero(bufred_d, &red_s, stream);
    }
    cuda_global_reduce_add(bufred, bufred_d, 1, stream);

    return bufred[0];
  }

   /**
   * Fortran wrapper glsum
   * Sum a vector of length n
   */
  real cuda_glsum(void *a, int *n, cudaStream_t stream) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    cuda_redbuf_check_alloc(nb);
    if ( *n > 0) {
      glsum_kernel<real>
        <<<nblcks, nthrds, 0, stream>>>((real *) a,
                                        (real *) bufred_d, *n);
      CUDA_CHECK(cudaGetLastError());
      reduce_kernel<real><<<1, 1024, 0, stream>>> ((real *) bufred_d, nb);
      CUDA_CHECK(cudaGetLastError());
    }
    else {
      cuda_rzero(bufred_d, &red_s, stream);
    }

    cuda_global_reduce_add(bufred, bufred_d, 1, stream);

    return bufred[0];
  }

  /**
   * Fortran wrapper absval
   * Take the abs value of a vector of length n
   */
  void cuda_absval(void *a, int *n, cudaStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    absval_kernel<real>
      <<<nblcks, nthrds,0, stream>>>((real *) a, * n);
    CUDA_CHECK(cudaGetLastError());

  }

  // ======================================================================== //
  // Point-wise operations.

  /** Fortran wrapper for pwmax_vec2
   *
   * Compute the maximum of two vectors \f$ a = \max(a, b) \f$
   */
  void cuda_pwmax_vec2(void *a, void *b, int *n, cudaStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      pwmax_vec2_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((real *)a, (real *)b, *n);
      CUDA_CHECK(cudaGetLastError());
  }

  /** Fortran wrapper for pwmax_vec3
   *
   * Compute the maximum of two vectors \f$ a = \max(b, c) \f$
   */
  void cuda_pwmax_vec3(void *a, void *b, void *c, int *n, cudaStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      pwmax_vec3_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((real *)a, (real *)b, (real *)c, *n);
      CUDA_CHECK(cudaGetLastError());
  }

  /** Fortran wrapper for pwmax_sca2
   *
   * Compute the maximum of vector and scalar \f$ a = \max(a, c) \f$
   */
  void cuda_pwmax_sca2(void *a, real *c, int *n, cudaStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      pwmax_sca2_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((real *)a, *c, *n);
      CUDA_CHECK(cudaGetLastError());
  }

  /** Fortran wrapper for pwmax_sca3
   *
   * Compute the maximum of vector and scalar \f$ a = \max(b, c) \f$
   */
  void cuda_pwmax_sca3(void *a, void *b, real *c, int *n, cudaStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      pwmax_sca3_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((real *)a, (real *)b, *c, *n);
      CUDA_CHECK(cudaGetLastError());
  }

  /** Fortran wrapper for pwmin_vec2
   *
   * Compute the minimum of two vectors \f$ a = \min(a, b) \f$
   */
  void cuda_pwmin_vec2(void *a, void *b, int *n, cudaStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      pwmin_vec2_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((real *)a, (real *)b, *n);
      CUDA_CHECK(cudaGetLastError());
  }

  /** Fortran wrapper for pwmin_vec3
   *
   * Compute the minimum of two vectors \f$ a = \min(b, c) \f$
   */
  void cuda_pwmin_vec3(void *a, void *b, void *c, int *n, cudaStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      pwmin_vec3_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((real *)a, (real *)b, (real *)c, *n);
      CUDA_CHECK(cudaGetLastError());
  }

  /** Fortran wrapper for pwmin_sca2
   *
   * Compute the minimum of vector and scalar \f$ a = \min(a, c) \f$
   */
  void cuda_pwmin_sca2(void *a, real *c, int *n, cudaStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      pwmin_sca2_kernel<real><<<nblcks, nthrds, 0, stream>>>((real *)a, *c, *n);
      CUDA_CHECK(cudaGetLastError());
  }

  /** Fortran wrapper for pwmin_sca3
   *
   * Compute the minimum of vector and scalar \f$ a = \min(b, c) \f$
   */
  void cuda_pwmin_sca3(void *a, void *b, real *c, int *n, cudaStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      pwmin_sca3_kernel<real><<<nblcks, nthrds, 0, stream>>>
        ((real *)a, (real *)b, *c, *n);
      CUDA_CHECK(cudaGetLastError());
  }

  // ======================================================================== //
  
  /** Fortran wrapper for iadd
   * Add a scalar to vector \f$ a_i = a_i + c \f$
   */
  void cuda_iadd(void *a, int *c, int *n, cudaStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    cadd_kernel<int><<<nblcks, nthrds, 0, stream>>>
      ((int *) a, *c, *n);
    CUDA_CHECK(cudaGetLastError());

  }

} /* extern "C" */
