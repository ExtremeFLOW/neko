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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "math_kernel.h"

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

#ifdef HAVE_RCCL
#include <math/bcknd/device/device_nccl_reduce.h>
#include <math/bcknd/device/device_nccl_op.h>
#endif

  /** Fortran wrapper for copy
   * Copy a vector \f$ a = b \f$
   */
  void hip_copy(void *a, void *b, int *n, hipStream_t strm) {
    HIP_CHECK(hipMemcpyAsync(a, b, (*n) * sizeof(real),
                             hipMemcpyDeviceToDevice, strm));
  }

  /** Fortran wrapper for masked copy
   * Copy a vector \f$ a(mask) = b(mask) \f$
   */
  void hip_masked_copy(void *a, void *b, void *mask,
                       int *n, int *m, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(masked_copy_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (int *) mask, *n, *m);

    HIP_CHECK(hipGetLastError());

  }

  /** Fortran wrapper for masked gather copy
   * Copy a vector \f$ a = b(mask) \f$
   */
  void hip_masked_gather_copy(void *a, void *b, void *mask,
                              int *n, int *m, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(masked_gather_copy_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (int *) mask, *n, *m);
    HIP_CHECK(hipGetLastError());

  }

  /** Fortran wrapper for masked scatter copy
   * Copy a vector \f$ a(mask) = b \f$
   */
  void hip_masked_scatter_copy(void *a, void *b, void *mask,
                               int *n, int *m, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(masked_scatter_copy_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (int *) mask, *n, *m);

    HIP_CHECK(hipGetLastError());

  }

  /** Fortran wrapper for masked atomic reduction
   * update a vector \f$ a += b(mask) \f$ where mask is not unique
   */
  void hip_masked_atomic_reduction(void *a, void *b, void *mask,
                                   int *n, int *m, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(masked_atomic_reduction_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (int *) mask, *n, *m);

    HIP_CHECK(hipGetLastError());

  }

  /** Fortran wrapper for cfill_mask
   * Fill a scalar to vector \f$ a_i = s, for i \in mask \f$
   */
  void hip_cfill_mask(void* a, real* c, int* size, void* mask, int* mask_size,
                      hipStream_t strm) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*mask_size) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cfill_mask_kernel<real>),
                       nblcks, nthrds, 0, strm, (real*)a,
                       *c, *size, (int*)mask, *mask_size);

    HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for rzero
   * Zero a real vector
   */
  void hip_rzero(void *a, int *n, hipStream_t strm) {
    HIP_CHECK(hipMemsetAsync(a, 0, (*n) * sizeof(real), strm));
  }

  /** Fortran wrapper for cmult
   * Multiplication by constant c \f$ a = c \cdot a \f$
   */
  void hip_cmult(void *a, real *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cmult_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a, *c, *n);
    HIP_CHECK(hipGetLastError());

  }

  /** Fortran wrapper for cmult
   * Multiplication by constant c \f$ a = c \cdot b \f$
   */
  void hip_cmult2(void *a, void *b, real *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cmult2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,(real *) b, *c, *n);
    HIP_CHECK(hipGetLastError());

  }

  /** Fortran wrapper for cdiv
   * Division of constant c by array \f$ a = c / a \f$
   */
  void hip_cdiv(void *a, real *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cdiv_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a, *c, *n);
    HIP_CHECK(hipGetLastError());

  }

  /** Fortran wrapper for cdiv
   * Division of constant c by \f$ a = c / b \f$
   */
  void hip_cdiv2(void *a, void *b, real *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cdiv2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,(real *) b, *c, *n);
    HIP_CHECK(hipGetLastError());

  }

  /** Fortran wrapper for cadd
   * Add a scalar to vector \f$ a = \sum a_i + s \f$
   */
  void hip_radd(void *a, real *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cadd_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a, *c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for cadd2
   * Add a scalar to a vector \f$ a = b + s \f$
   */
  void hip_cadd2(void *a, void *b, real *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cadd2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a, (real *) b, *c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for cfill
   * Multiplication by constant c \f$ a = c \cdot a \f$
   */
  void hip_cfill(void *a, real *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    if (*n > 0) {
      hipLaunchKernelGGL(HIP_KERNEL_NAME(cfill_kernel<real>),
                         nblcks, nthrds, 0, strm, (real *) a, *c, *n);
      HIP_CHECK(hipGetLastError());
    }

  }

  /**
   * Fortran wrapper for add2
   * Vector addition \f$ a = a + b \f$
   */
  void hip_add2(void *a, void *b, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a, (real *) b, *n);
    HIP_CHECK(hipGetLastError());

  }

  /**
   * Fortran wrapper for add4
   * Vector addition \f$ a = b + c + d\f$
   */
  void hip_add3(void *a, void *b, void *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add3_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add4
   * Vector addition \f$ a = b + c + d\f$
   */
  void hip_add4(void *a, void *b, void *c, void *d, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add4_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, (real *) d, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add2s1
   * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
   * (multiplication on first argument)
   */
  void hip_add2s1(void *a, void *b, real *c1, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2s1_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, *c1, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add2s2
   * Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
   * (multiplication on second argument)
   */
  void hip_add2s2(void *a, void *b, real *c1, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2s2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a, (real *) b,
                       *c1, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add2s2
   * Vector addition with scalar multiplication
   * \f$ x = x + c_1 p1 + c_2p2 + ... + c_jpj \f$
   * (multiplication on second argument)
   */
  void hip_add2s2_many(void *x, void **p, void *alpha, int *j, int *n,
                       hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2s2_many_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) x,
                       (const real **) p, (real *) alpha, *j, *n);
    HIP_CHECK(hipGetLastError());

  }

  /**
   * Fortran wrapper for addsqr2s2
   * Vector addition with scalar multiplication \f$ a = a + c_1 (b * b) \f$
   * (multiplication on second argument)
   */
  void hip_addsqr2s2(void *a, void *b, real *c1, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(addsqr2s2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, *c1, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add3s2
   * Vector addition with scalar multiplication \f$ a = c_1 b + c_2 c \f$
   * (multiplication on second argument)
   */
  void hip_add3s2(void *a, void *b, void *c, real *c1, real *c2, int *n,
                  hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add3s2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, *c1, *c2, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add4s3
   * Vector addition with scalar multiplication \f$ a = c_1 b + c_2 c + c_3 d \f$
   * (multiplication on second argument)
   */
  void hip_add4s3(void *a, void *b, void *c, void *d, real *c1, real *c2,
                  real *c3, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add4s3_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, (real *) d, *c1, *c2, *c3, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add5s4
   * Vector addition with scalar multiplication \f$ a = a + c_1 b + c_2 c + c_3 d + c_4 e \f$
   * (multiplication on second argument)
   */
  void hip_add5s4(void *a, void *b, void *c, void *d, void *e, real *c1,
                  real *c2, real *c3, real *c4, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add5s4_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, (real *) d, (real *) e,
                       *c1, *c2, *c3, *c4, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for invcol1
   * Invert a vector \f$ a = 1 / a \f$
   */
  void hip_invcol1(void *a, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(invcol1_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for invcol2
   * Vector division \f$ a = a / b \f$
   */
  void hip_invcol2(void *a, void *b, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(invcol2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a, (real *) b, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for invcol3
   * Vector division \f$ a = b / c \f$
   */
  void hip_invcol3(void *a, void *b, void *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(invcol3_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b,  (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for col2
   * Vector multiplication with 2 vectors \f$ a = a \cdot b \f$
   */
  void hip_col2(void *a, void *b, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(col2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for col3
   * Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
   */
  void hip_col3(void *a, void *b, void *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(col3_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for subcol3
   * Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
   */
  void hip_subcol3(void *a, void *b, void *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(subcol3_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for sub2
   * Vector subtraction \f$ a = a - b \f$
   */
  void hip_sub2(void *a, void *b, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sub2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for sub3
   * Vector subtraction \f$ a = b - c \f$
   */
  void hip_sub3(void *a, void *b, void *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sub3_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for addcol3
   * \f$ a = a + b * c \f$
   */
  void hip_addcol3(void *a, void *b, void *c, int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(addcol3_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for addcol4
   * \f$ a = a + b * c * d \f$
   */
  void hip_addcol4(void *a, void *b, void *c, void *d, int *n,
                   hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(addcol4_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, (real *) d, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for addcol3s2
   * \f$ a = a + s(b * c) \f$
   */
  void hip_addcol3s2(void *a, void *b, void *c, real *s, int *n,
                     hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(addcol3s2_kernel<real>),
                       nblcks, nthrds, 0, strm, (real *) a,
                       (real *) b, (real *) c, *s, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for vdot3
   * \f$ dot = u \cdot v \f$
   */
  void hip_vdot3(void *dot, void *u1, void *u2, void *u3,
                 void *v1, void *v2, void *v3, int *n,
                 hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(vdot3_kernel<real>),
                       nblcks, nthrds, 0, strm,
                       (real *) dot, (real *) u1, (real *) u2, (real *) u3,
                       (real *) v1, (real *) v2, (real *) v3, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for vdot3
   * \f$ u = v \times w \f$
   */
  void hip_vcross(void *u1, void *u2, void *u3,
                  void *v1, void *v2, void *v3,
                  void *w1, void *w2, void *w3,
                  int *n, hipStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(vcross_kernel<real>),
                       nblcks, nthrds, 0, strm,
                       (real *) u1, (real *) u2, (real *) u3,
                       (real *) v1, (real *) v2, (real *) v3,
                       (real *) w1, (real *) w2, (real *) w3, *n);
    HIP_CHECK(hipGetLastError());
  }


  /*
   * Reduction buffer
   */
  int red_s = 0;
  real * bufred = NULL;
  real * bufred_d = NULL;

  void hip_redbuf_check_alloc(int nb) {
    if ( nb >= red_s) {
      red_s = nb+1;
      if (bufred != NULL) {
        HIP_CHECK(hipHostFree(bufred));
        HIP_CHECK(hipFree(bufred_d));
      }
      HIP_CHECK(hipHostMalloc(&bufred,red_s*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&bufred_d, red_s*sizeof(real)));
    }
  }

  /**
   * Global additive reduction
   */
  void hip_global_reduce_add(real * bufred, void * bufred_d, int n,
                             const hipStream_t stream) {
    #ifdef HAVE_RCCL
        device_nccl_allreduce(bufred_d, bufred_d, n, sizeof(real),
                              DEVICE_NCCL_SUM, stream);
        HIP_CHECK(hipMemcpyAsync(bufred, bufred_d, n*sizeof(real),
                                 hipMemcpyDeviceToHost, stream));
        HIP_CHECK(hipStreamSynchronize(stream));
    #elif HAVE_DEVICE_MPI
        HIP_CHECK(hipStreamSynchronize(stream));
        device_mpi_allreduce(bufred_d, bufred, n, sizeof(real), DEVICE_MPI_SUM);
    #else
        HIP_CHECK(hipMemcpyAsync(bufred, bufred_d, n*sizeof(real),
                                 hipMemcpyDeviceToHost, stream));
        HIP_CHECK(hipStreamSynchronize(stream));
    #endif
  }

  /**
   * Fortran wrapper vlsc3
   * Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
   */
  real hip_vlsc3(void *u, void *v, void *w, int *n, hipStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    hip_redbuf_check_alloc(nb);


    hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) u, (real *) v,
                       (real *) w, bufred_d, *n);
    HIP_CHECK(hipGetLastError());
    hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0, stream, bufred_d, nb);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpyAsync(bufred, bufred_d, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));

    return bufred[0];
  }

  /**
   * Fortran wrapper glsc3
   * Weighted inner product \f$ a^T b c \f$
   */
  real hip_glsc3(void *a, void *b, void *c, int *n, hipStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    hip_redbuf_check_alloc(nb);

    if (*n > 0) {
      hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_kernel<real>),
                         nblcks, nthrds, 0, stream,
                         (real *) a, (real *) b,
                         (real *) c, bufred_d, *n);
      HIP_CHECK(hipGetLastError());
      hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                         1, 1024, 0, stream, bufred_d, nb);
      HIP_CHECK(hipGetLastError());
    }
    else {
      hip_rzero(bufred_d, &red_s, stream);
    }
    hip_global_reduce_add(bufred, bufred_d, 1, stream);

    return bufred[0];
  }

  /**
   * Fortran wrapper for doing an reduction to an array
   * Weighted inner product \f$ w^T v(n,1:j) c \f$
   */
  void hip_glsc3_many(real *h, void * w, void *v,void *mult, int *j, int *n,
                      hipStream_t stream){
    int pow2 = 1;
    while(pow2 < (*j)){
      pow2 = 2*pow2;
    }
    const int nt = 1024/pow2;
    const dim3 nthrds(pow2, nt, 1);
    const dim3 nblcks(((*n)+nt - 1)/nt, 1, 1);
    const dim3 nthrds_red(1024,1,1);
    const dim3 nblcks_red( (*j),1,1);
    const int nb = ((*n) + nt - 1)/nt;

    hip_redbuf_check_alloc((*j)*nb);

    if (*n > 0) {
      hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_many_kernel<real>),
                         nblcks, nthrds, 0, stream,
                         (const real *) w, (const real **) v,
                         (const real *)mult, bufred_d, *j, *n);
      HIP_CHECK(hipGetLastError());

      hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_reduce_kernel<real>),
                         nblcks_red, nthrds_red, 0, stream,
                         bufred_d, nb, *j);
      HIP_CHECK(hipGetLastError());
    }
    else {
      hip_rzero(bufred_d, &red_s, stream);
    }
    hip_global_reduce_add(h, bufred_d, (*j), stream);
  }

  /**
   * Fortran wrapper glsc2
   * Weighted inner product \f$ a^T b \f$
   */
  real hip_glsc2(void *a, void *b, int *n, hipStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;


    hip_redbuf_check_alloc(nb);

    if( *n > 0) {
      hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc2_kernel<real>),
                         nblcks, nthrds, 0, stream,
                         (real *) a, (real *) b, bufred_d, *n);
      HIP_CHECK(hipGetLastError());
      hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                         1, 1024, 0, stream, bufred_d, nb);
      HIP_CHECK(hipGetLastError());
    }
    else {
      hip_rzero(bufred_d, &red_s, stream);
    }
    hip_global_reduce_add(bufred, bufred_d, 1, stream);

    return bufred[0];
  }

  /**
   * Fortran wrapper glsubnorm2
   * Weighted inner product \f$ a^T b \f$
   */
  real hip_glsubnorm2(void* a, void* b, int* n, hipStream_t stream) {

      const dim3        nthrds(1024, 1, 1);
      const dim3        nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
      const int         nb     = ((*n) + 1024 - 1) / 1024;

      hip_redbuf_check_alloc(nb);

      if (*n > 0) {
        hipLaunchKernelGGL(HIP_KERNEL_NAME(glsubnorm2_kernel<real>),
                           nblcks, nthrds, 0, stream,
                           (real*)a, (real*)b, bufred_d, *n);
        HIP_CHECK(hipGetLastError());
        hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                           1, 1024, 0, stream, bufred_d, nb);
        HIP_CHECK(hipGetLastError());
      }
      else {
        hip_rzero(bufred_d, &red_s, stream);
      }
      hip_global_reduce_add(bufred, bufred_d, 1, stream);

      return bufred[0];
  }

  /**
   * Fortran wrapper glsum
   * Sum a vector of length n
   */
  real hip_glsum(void *a, int *n, hipStream_t stream) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    hip_redbuf_check_alloc(nb);
    if( *n > 0) {
      hipLaunchKernelGGL(HIP_KERNEL_NAME(glsum_kernel<real>),
                         nblcks, nthrds, 0, stream,
                         (real *) a, bufred_d, *n);
      HIP_CHECK(hipGetLastError());
      hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                         1, 1024, 0, stream, bufred_d, nb);
      HIP_CHECK(hipGetLastError());
    }
    else {
      hip_rzero(bufred_d, &red_s, stream);
    }

    hip_global_reduce_add(bufred, bufred_d, 1, stream);

    return bufred[0];
  }

  /**
   * Fortran wrapper for absval
   * \f$ a = abs(a) \f$
   */
  void hip_absval(void *a, int *n, hipStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(absval_kernel<real>),
                       nblcks, nthrds, 0, stream, (real *) a, *n);
    HIP_CHECK(hipGetLastError());

}

  // ======================================================================== //
  // Point-wise operations.

  /** Fortran wrapper for pwmax_vec2
   *
   * Compute the maximum of two vectors \f$ a = \max(a, b) \f$
   */
  void hip_pwmax_vec2(void* a, void* b, int* n, hipStream_t stream) {

    const dim3        nthrds(1024, 1, 1);
    const dim3        nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(pwmax_vec2_kernel<real>),
                       nblcks, nthrds, 0, stream, (real*)a, (real*)b, *n);
    HIP_CHECK(hipGetLastError());
}

  /** Fortran wrapper for pwmax_vec3
   *
   * Compute the maximum of two vectors \f$ a = \max(b, c) \f$
   */
  void hip_pwmax_vec3(void *a, void *b, void *c, int *n, hipStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(pwmax_vec3_kernel<real>),
                         nblcks, nthrds, 0, stream,
                         (real *)a, (real *)b, (real *)c, *n);
      HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for pwmax_sca2
   *
   * Compute the maximum of vector and scalar \f$ a = \max(a, c) \f$
   */
  void hip_pwmax_sca2(void *a, real *c, int *n, hipStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(pwmax_sca2_kernel<real>),
                         nblcks, nthrds, 0, stream, (real *)a, *c, *n);
      HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for pwmax_sca3
   *
   * Compute the maximum of vector and scalar \f$ a = \max(b, c) \f$
   */
  void hip_pwmax_sca3(void *a, void *b, real *c, int *n, hipStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(pwmax_sca3_kernel<real>),
                         nblcks, nthrds, 0, stream,
                         (real *)a, (real *)b, *c, *n);
      HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for pwmin_vec2
   *
   * Compute the minimum of two vectors \f$ a = \min(a, b) \f$
   */
  void hip_pwmin_vec2(void *a, void *b, int *n, hipStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(pwmin_vec2_kernel<real>),
                         nblcks, nthrds, 0, stream,
                         (real *)a, (real *)b, *n);
      HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for pwmin_vec3
   *
   * Compute the minimum of two vectors \f$ a = \min(b, c) \f$
   */
  void hip_pwmin_vec3(void *a, void *b, void *c, int *n, hipStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(pwmin_vec3_kernel<real>),
                         nblcks, nthrds, 0, stream,
                        (real *)a, (real *)b, (real *)c, *n);
      HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for pwmin_sca2
   *
   * Compute the minimum of vector and scalar \f$ a = \min(a, c) \f$
   */
  void hip_pwmin_sca2(void *a, real *c, int *n, hipStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(pwmin_sca2_kernel<real>),
                         nblcks, nthrds, 0, stream,
                        (real *)a, *c, *n);
      HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for pwmin_sca3
   *
   * Compute the minimum of vector and scalar \f$ a = \min(b, c) \f$
   */
  void hip_pwmin_sca3(void *a, void *b, real *c, int *n, hipStream_t stream) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(pwmin_sca3_kernel<real>),
                         nblcks, nthrds, 0, stream,
                         (real *)a, (real *)b, *c, *n);
      HIP_CHECK(hipGetLastError());
  }

  // ======================================================================== //

  /** Fortran wrapper for cadd
   * Add a scalar to vector \f$ a = \sum a_i + s \f$
   */
  void hip_iadd(void *a, int *c, int *n, hipStream_t stream) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cadd_kernel<int>),
                       nblcks, nthrds, 0, stream,
                       (int *) a, *c, *n);
    HIP_CHECK(hipGetLastError());
  }

} /* extern "C" */
