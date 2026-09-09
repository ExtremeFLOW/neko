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

#include <stdio.h>
#include <hip/hip_runtime.h>
#include "bicgstab_kernel.h"
#include <device/device_config.h>
#include <device/hip/check.h>
#include <device/hip/buffer.h>

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

#ifdef HAVE_RCCL
#include <math/bcknd/device/device_nccl_reduce.h>
#include <math/bcknd/device/device_nccl_op.h>
#endif

  /**
   * Reduction buffer, owned by the device layer and released on
   * device teardown (hip_buffer_free_all in hip_finalize)
   */
  hip_buffer_t bicgstab_redbuf = HIP_BUFFER_INIT;

  /**
   * Sum @a count contiguous values across all ranks and stage the result
   * in the pinned host buffer
   */
  static void bicgstab_reduce(real_xp *host, real_xp *dev, const int count,
                              hipStream_t stream) {
#ifdef HAVE_RCCL
    device_nccl_allreduce(dev, dev, count, sizeof(real_xp),
                          DEVICE_NCCL_SUM, stream);
    HIP_CHECK(hipMemcpyAsync(host, dev, count * sizeof(real_xp),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#elif HAVE_DEVICE_MPI
    HIP_CHECK(hipStreamSynchronize(stream));
    device_mpi_allreduce(dev, host, count, sizeof(real_xp), DEVICE_MPI_SUM);
#else
    HIP_CHECK(hipMemcpyAsync(host, dev, count * sizeof(real_xp),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#endif
  }

  /**
   * Fortran wrapper for the BiCGStab search direction update
   * \f$ p = r + \beta (p - \omega v) \f$
   */
  void hip_bicgstab_update_p(void *p, void *r, void *v, real *beta,
                             real *omega, int *n) {

    if (*n <= 0)
      return;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hipLaunchKernelGGL(HIP_KERNEL_NAME(bicgstab_update_p_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) p, (real *) r, (real *) v,
                       *beta, *omega, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for a weighted inner product and squared norm,
   * \f$ (a^T M b, b^T M b) \f$
   */
  void hip_bicgstab_product_and_norm(void *a, void *b, void *mult,
                                     real_xp *res, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1) / 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hip_buffer_reserve(&bicgstab_redbuf,
                       2 * (nb > 1 ? nb : 1) * sizeof(real_xp));
    real_xp *host = (real_xp *) bicgstab_redbuf.host;
    real_xp *dev = (real_xp *) bicgstab_redbuf.dev;

    if (*n > 0) {
      hipLaunchKernelGGL(
        HIP_KERNEL_NAME(bicgstab_product_and_norm_kernel<real, real_xp>),
        nblcks, nthrds, 0, stream,
        (real *) a, (real *) b, (real *) mult, dev, *n);
      HIP_CHECK(hipGetLastError());
      hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_reduce_kernel<real_xp>),
                         2, 1024, 0, stream, dev, nb, 2);
      HIP_CHECK(hipGetLastError());
    }
    else {
      HIP_CHECK(hipMemsetAsync(dev, 0, 2 * sizeof(real_xp), stream));
    }

    bicgstab_reduce(host, dev, 2, stream);
    res[0] = host[0];
    res[1] = host[1];
  }

  /**
   * Fortran wrapper for BiCGStab part 1, \f$ s = r - \alpha v \f$,
   * returning \f$ s^T M s \f$
   */
  real_xp hip_bicgstab_part1(void *s, void *r, void *v, void *mult,
                             real *alpha, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1) / 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hip_buffer_reserve(&bicgstab_redbuf,
                       2 * (nb > 1 ? nb : 1) * sizeof(real_xp));
    real_xp *host = (real_xp *) bicgstab_redbuf.host;
    real_xp *dev = (real_xp *) bicgstab_redbuf.dev;

    if (*n > 0) {
      hipLaunchKernelGGL(HIP_KERNEL_NAME(bicgstab_part1_kernel<real, real_xp>),
                         nblcks, nthrds, 0, stream,
                         (real *) s, (real *) r, (real *) v,
                         (real *) mult, dev, *alpha, *n);
      HIP_CHECK(hipGetLastError());
      hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real_xp>),
                         1, 1024, 0, stream, dev, nb);
      HIP_CHECK(hipGetLastError());
    }
    else {
      HIP_CHECK(hipMemsetAsync(dev, 0, sizeof(real_xp), stream));
    }

    bicgstab_reduce(host, dev, 1, stream);
    return host[0];
  }

  /**
   * Fortran wrapper for BiCGStab part 2,
   * \f$ x = x + \alpha \hat{p} + \omega \hat{s} \f$ and
   * \f$ r = s - \omega t \f$, returning
   * \f$ (r^T M r, f^T M r) \f$
   */
  void hip_bicgstab_part2(void *x, void *r, void *p_hat, void *s_hat,
                          void *s, void *t, void *f, void *mult,
                          real *alpha, real *omega, real_xp *res, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1) / 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hip_buffer_reserve(&bicgstab_redbuf,
                       2 * (nb > 1 ? nb : 1) * sizeof(real_xp));
    real_xp *host = (real_xp *) bicgstab_redbuf.host;
    real_xp *dev = (real_xp *) bicgstab_redbuf.dev;

    if (*n > 0) {
      hipLaunchKernelGGL(HIP_KERNEL_NAME(bicgstab_part2_kernel<real, real_xp>),
                         nblcks, nthrds, 0, stream,
                         (real *) x, (real *) r,
                         (real *) p_hat, (real *) s_hat,
                         (real *) s, (real *) t, (real *) f,
                         (real *) mult, dev, *alpha, *omega, *n);
      HIP_CHECK(hipGetLastError());
      hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_reduce_kernel<real_xp>),
                         2, 1024, 0, stream, dev, nb, 2);
      HIP_CHECK(hipGetLastError());
    }
    else {
      HIP_CHECK(hipMemsetAsync(dev, 0, 2 * sizeof(real_xp), stream));
    }

    bicgstab_reduce(host, dev, 2, stream);
    res[0] = host[0];
    res[1] = host[1];
  }
}
