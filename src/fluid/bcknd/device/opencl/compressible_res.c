/*
 Copyright (c) 2025, The Neko Authors
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

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <device/device_config.h>
#include <device/opencl/jit.h>
#include <device/opencl/prgm_lib.h>
#include <device/opencl/check.h>
#include <common/neko_log.h>

#include "compressible_res_kernel.cl.h"

/*
 * Tuning candidates of the fused divergence
 *
 * The kernel's geometry is measured rather than decided: the tuner below
 * times the candidates on the first call at each order and keeps the
 * fastest. Candidate c selects, in its bits, which of the point state
 * (bit 0) and the geometric factors (bit 1) of a work-item's points are held
 * in registers across the passes rather than re-read from global memory,
 * and whether three fluxes or one are staged in local memory per pass
 * (bit 2). The CUDA and HIP kernels also offer a narrower work-group in
 * bit 3; here the work-group is capped at 256 work-items already, so that
 * bit is ignored and the candidate takes the geometry of the one without it.
 * Candidates that round to the same geometry name the same kernel and the
 * tuner skips the repeats. The kernel names are generated from these rules,
 * see the instantiation list in compressible_res_kernel.cl.
 */
#define INVISCID_DIV_CANDIDATES 16
#define INVISCID_DIV_WARP 64
#define INVISCID_DIV_MAX_BLOCK 256
#define INVISCID_DIV_MAX_SMEM 49152
#define INVISCID_DIV_MAX_REG_PTS 4
#define INVISCID_DIV_TUNE_WARMUP 20

/* Work-group size: the points of one chunk, as few chunks of at most
   INVISCID_DIV_MAX_BLOCK work-items as needed, rounded up to whole waves */
static int inviscid_div_chunks(const int lx) {
  const int npts = lx * lx * lx;
  const int nchunks = (npts + INVISCID_DIV_MAX_BLOCK - 1)
    / INVISCID_DIV_MAX_BLOCK;
  const int per_chunk = (npts + nchunks - 1) / nchunks;
  return ((per_chunk + INVISCID_DIV_WARP - 1) / INVISCID_DIV_WARP)
    * INVISCID_DIV_WARP;
}

/* Points per work-item */
static int inviscid_div_pts(const int lx) {
  return (lx * lx * lx + inviscid_div_chunks(lx) - 1) / inviscid_div_chunks(lx);
}

/* Fluxes staged per pass. The kernels are always compiled in double
   precision (see scripts/cltostring.sh), whatever the host's real is */
static int inviscid_div_stage(const int lx, const int c) {
  return (((c >> 2) & 1) == 0 &&
          (3 * lx * lx * lx + 3 * lx * lx) * (int) sizeof(double)
          <= INVISCID_DIV_MAX_SMEM) ? 3 : 1;
}

/* Register choice: bit 0 the state, bit 1 the geometric factors, neither
   offered beyond INVISCID_DIV_MAX_REG_PTS points per work-item */
static int inviscid_div_reg(const int lx, const int c) {
  return (inviscid_div_pts(lx) <= INVISCID_DIV_MAX_REG_PTS) ? (c & 3) : 0;
}

/* Does candidate c name the same kernel as a lower candidate? */
static int inviscid_div_cand_dup(const int lx, const int c) {
  for (int d = 0; d < c; d++) {
    if (inviscid_div_stage(lx, d) == inviscid_div_stage(lx, c) &&
        inviscid_div_reg(lx, d) == inviscid_div_reg(lx, c)) {
      return 1;
    }
  }
  return 0;
}

/* Describe candidate c for the log */
static const char *inviscid_div_cand_desc(const int lx, const int c,
                                          char *buf, const int len) {
  static const char *regs[] = { "none", "state", "metrics", "both" };
  snprintf(buf, len, "%d wi, %d/pass, regs: %s", inviscid_div_chunks(lx),
           inviscid_div_stage(lx, c), regs[inviscid_div_reg(lx, c)]);
  return buf;
}

/* Candidate pinned by NEKO_INVISCID_DIV, or -1 to sweep */
static int neko_inviscid_div_pin(void) {
  const char *v = getenv("NEKO_INVISCID_DIV");
  int c;

  if (v == NULL) {
    return -1;
  }

  c = atoi(v);
  if (c < 0 || c >= INVISCID_DIV_CANDIDATES) {
    c = 0;
  }
  return c;
}

/* Sampling of the sweep, as for the CUDA and HIP operator tuners */
static int neko_tune_rounds(void) {
  const char *v = getenv("NEKO_TUNE_ROUNDS");
  int n = (v != NULL) ? atoi(v) : 3;

  return (n < 1) ? 1 : n;
}

static int neko_tune_iters(void) {
  const char *v = getenv("NEKO_TUNE_ITERS");
  int n = (v != NULL) ? atoi(v) : 100;

  return (n < 1) ? 1 : n;
}

/* Enqueue candidate c at order lx on QUEUE, one work-group per element */
static void inviscid_div_launch(const int lx, const int c, void **bufs,
                                const int nel, cl_command_queue queue,
                                cl_event *event) {
  cl_int err;
  char name[96];

  snprintf(name, sizeof(name),
           "compressible_res_inviscid_div_kernel_lx%d_ch%d_nc%d_r%d",
           lx, inviscid_div_chunks(lx), inviscid_div_stage(lx, c),
           inviscid_div_reg(lx, c));
  cl_kernel kernel = clCreateKernel(compressible_res_program, name, &err);
  CL_CHECK(err);

  for (int i = 0; i < 27; i++) {
    CL_CHECK(clSetKernelArg(kernel, i, sizeof(cl_mem), (void *) &bufs[i]));
  }

  const size_t local_item_size = inviscid_div_chunks(lx);
  const size_t global_item_size = local_item_size * nel;

  CL_CHECK(clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_item_size,
                                  &local_item_size, 0, NULL, event));
  CL_CHECK(clReleaseKernel(kernel));
}

/*
 * Time every distinct candidate at this order on the profiling queue and
 * return the fastest: interleaved rounds, min reduced, after a warm-up of
 * every candidate, as the CUDA and HIP operator tuners sample. The
 * sampling and the hand-over between the queues follow opencl_dudxyz().
 * NEKO_INVISCID_DIV pins a candidate, NEKO_TUNE_ROUNDS and NEKO_TUNE_ITERS
 * set the sampling.
 *
 * Every candidate computes the same sums in the same order, so the tuner
 * leaves the output of the chosen one in place merely by convention: it
 * stands in for the evaluation it was called for.
 */
static int tune_inviscid_div(const int lx, void **bufs, const int nel) {
  char neko_log_buf[80];
  char desc[48];
  double time[INVISCID_DIV_CANDIDATES];
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int pin = neko_inviscid_div_pin();
  int best = 0;

  sprintf(neko_log_buf, "Autotune inviscid div (lx: %d)", lx);
  log_section(neko_log_buf);

  if (pin >= 0) {
    inviscid_div_launch(lx, pin, bufs, nel, (cl_command_queue) glb_cmd_queue,
                        NULL);
    sprintf(neko_log_buf, "Set by env   : %d (%s)", pin,
            inviscid_div_cand_desc(lx, pin, desc, sizeof(desc)));
    log_message(neko_log_buf);
    log_end_section();
    return pin;
  }

  for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
    time[c] = HUGE_VAL;
  }

  cl_event perf_event, sync_event;
  cl_ulong start, end;

  CL_CHECK(clFinish(glb_cmd_queue));
  CL_CHECK(clEnqueueMarker(glb_cmd_queue, &sync_event));
  CL_CHECK(clEnqueueBarrier(prf_cmd_queue));
  CL_CHECK(clEnqueueWaitForEvents(prf_cmd_queue, 1, &sync_event));

  /* Warm every candidate before timing anything. Candidates that name the
     kernel of a lower one are skipped throughout */
  for (int i = 0; i < INVISCID_DIV_TUNE_WARMUP; i++) {
    for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
      if (!inviscid_div_cand_dup(lx, c)) {
        inviscid_div_launch(lx, c, bufs, nel, (cl_command_queue) prf_cmd_queue,
                            NULL);
      }
    }
  }
  CL_CHECK(clFinish(prf_cmd_queue));

  /* Interleaved rounds, best summed kernel time per candidate */
  for (int r = 0; r < rounds; r++) {
    for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
      if (inviscid_div_cand_dup(lx, c)) {
        continue;
      }
      double elapsed = 0.0;
      for (int i = 0; i < iters; i++) {
        inviscid_div_launch(lx, c, bufs, nel, (cl_command_queue) prf_cmd_queue,
                            &perf_event);
        CL_CHECK(clWaitForEvents(1, &perf_event));
        CL_CHECK(clGetEventProfilingInfo(perf_event,
                                         CL_PROFILING_COMMAND_START,
                                         sizeof(cl_ulong), &start, NULL));
        CL_CHECK(clGetEventProfilingInfo(perf_event,
                                         CL_PROFILING_COMMAND_END,
                                         sizeof(cl_ulong), &end, NULL));
        CL_CHECK(clReleaseEvent(perf_event));
        elapsed += (end - start) * 1.0e-6;
      }
      if (elapsed < time[c]) {
        time[c] = elapsed;
      }
    }
  }

  for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
    if (time[c] >= HUGE_VAL) {
      continue;
    }
    sprintf(neko_log_buf, "cand %2d (%s): %9.2f us/call", c,
            inviscid_div_cand_desc(lx, c, desc, sizeof(desc)),
            time[c] * 1000.0 / (double) iters);
    log_message(neko_log_buf);
    if (time[c] < time[best]) {
      best = c;
    }
  }

  inviscid_div_launch(lx, best, bufs, nel, (cl_command_queue) prf_cmd_queue,
                      NULL);
  CL_CHECK(clFinish(prf_cmd_queue));
  CL_CHECK(clEnqueueMarker(prf_cmd_queue, &sync_event));
  CL_CHECK(clEnqueueBarrier(glb_cmd_queue));
  CL_CHECK(clEnqueueWaitForEvents(glb_cmd_queue, 1, &sync_event));

  sprintf(neko_log_buf, "Chose        : %d (%s)", best,
          inviscid_div_cand_desc(lx, best, desc, sizeof(desc)));
  log_message(neko_log_buf);
  log_end_section();
  return best;
}

void compressible_res_part_visc_opencl(void *rhs_rho, void *rhs_m_x,
                                void *rhs_m_y, void *rhs_m_z, void *rhs_E,
                                void *visc_rho, void *visc_m_x, void *visc_m_y,
                                void *visc_m_z, void *visc_E,
                                void *mult, void *Binv, void *h1, int *n) {
  cl_int err;
  
  if (compressible_res_program == NULL)
    opencl_kernel_jit(compressible_res_kernel, (cl_program *) &compressible_res_program);

  cl_kernel kernel = clCreateKernel(compressible_res_program,
                                    "compressible_res_part_visc_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &rhs_rho));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &rhs_m_x));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &rhs_m_y));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &rhs_m_z));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &rhs_E));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &visc_rho));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &visc_m_x));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &visc_m_y));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &visc_m_z));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &visc_E));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &mult));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &Binv));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &h1));
  CL_CHECK(clSetKernelArg(kernel, 13, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
  CL_CHECK(clReleaseKernel(kernel));
}

/**
 * Work-group size of the fused divergence: the lx^3 points of an element
 * are split over as few chunks of at most 256 work-items as needed, rounded
 * up to 64. Must match the instantiation table in the kernel source.
 */
/**
 * Fused inviscid divergence of the compressible residual: one work-group
 * per element, see compressible_res_inviscid_div_kernel. The geometry is
 * autotuned on the first call at each order.
 */
void compressible_res_inviscid_div_opencl(void *rhs_rho, void *rhs_m_x,
                                          void *rhs_m_y, void *rhs_m_z,
                                          void *rhs_E, void *rho, void *m_x,
                                          void *m_y, void *m_z, void *E,
                                          void *p, void *u, void *v, void *w,
                                          void *drdx, void *dsdx, void *dtdx,
                                          void *drdy, void *dsdy, void *dtdy,
                                          void *drdz, void *dsdz, void *dtdz,
                                          void *dx, void *dy, void *dz,
                                          void *jacinv, int *nel, int *lx) {
  /* Candidate chosen by the tuner plus one, zero until it has run */
  static int autotune[17] = { 0 };
  void *bufs[27] = { rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E,
                     rho, m_x, m_y, m_z, E, p, u, v, w,
                     drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz,
                     dx, dy, dz, jacinv };

  if ((*lx) < 2 || (*lx) > 16) {
    fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
    exit(1);
  }

  if (compressible_res_program == NULL)
    opencl_kernel_jit(compressible_res_kernel, (cl_program *) &compressible_res_program);

  if (autotune[*lx] == 0) {
    autotune[*lx] = 1 + tune_inviscid_div(*lx, bufs, *nel);
  } else {
    inviscid_div_launch(*lx, autotune[*lx] - 1, bufs, *nel,
                        (cl_command_queue) glb_cmd_queue, NULL);
  }
}

void compressible_res_part_rk_sum_opencl(void *rho, void *m_x, void *m_y, void *m_z,
                                  void *E, void *k_rho_i, void *k_m_x_i,
                                  void *k_m_y_i, void *k_m_z_i, void *k_E_i,
                                  real *dt, real *c, int *n) {
  cl_int err;
  
  if (compressible_res_program == NULL)
    opencl_kernel_jit(compressible_res_kernel, (cl_program *) &compressible_res_program);
  
  cl_kernel kernel = clCreateKernel(compressible_res_program,
                                    "compressible_res_part_rk_sum_kernel", &err);
  CL_CHECK(err);

  CL_CHECK(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &rho));
  CL_CHECK(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &m_x));
  CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &m_y));
  CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &m_z));
  CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &E));
  CL_CHECK(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &k_rho_i));
  CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &k_m_x_i));
  CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &k_m_y_i));
  CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &k_m_z_i));
  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &k_E_i));
  CL_CHECK(clSetKernelArg(kernel, 10, sizeof(real), dt));
  CL_CHECK(clSetKernelArg(kernel, 11, sizeof(real), c));
  CL_CHECK(clSetKernelArg(kernel, 12, sizeof(int), n));
  
  const int nb = ((*n) + 256 - 1) / 256;
  const size_t global_item_size = 256 * nb;
  const size_t local_item_size = 256;

  CL_CHECK(clEnqueueNDRangeKernel((cl_command_queue) glb_cmd_queue, kernel, 1,
                                  NULL, &global_item_size, &local_item_size,
                                  0, NULL, NULL));  
  CL_CHECK(clReleaseKernel(kernel));
}
