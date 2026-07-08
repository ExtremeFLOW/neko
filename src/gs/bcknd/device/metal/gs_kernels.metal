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
 * Metal Shading Language kernels for gather-scatter operations.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Emulated float atomic min / max via compare-and-swap
 * (Metal only supports atomic_fetch_min/max for integer types)
 */

inline void atomic_float_min(device atomic_float *addr, float val) {
  float old = atomic_load_explicit(addr, memory_order_relaxed);
  while (val < old) {
    if (atomic_compare_exchange_weak_explicit(
            addr, &old, val,
            memory_order_relaxed, memory_order_relaxed))
      break;
  }
}

inline void atomic_float_max(device atomic_float *addr, float val) {
  float old = atomic_load_explicit(addr, memory_order_relaxed);
  while (val > old) {
    if (atomic_compare_exchange_weak_explicit(
            addr, &old, val,
            memory_order_relaxed, memory_order_relaxed))
      break;
  }
}

/*
 * Gather kernels
 */

/**
 * Gather kernel — addition
 * v(dg(i)) = sum of u(gd(k..k+blk_len-1))  for block i
 * Facet pairs: u(gd(i)) + u(gd(i+1))  or simple copy
 */
kernel void gather_kernel_add(device float *v[[ buffer(0) ]],
                              constant int &m[[ buffer(1) ]],
                              constant int &o[[ buffer(2) ]],
                              device const int *dg[[ buffer(3) ]],
                              device const float *u[[ buffer(4) ]],
                              constant int &n[[ buffer(5) ]],
                              device const int *gd[[ buffer(6) ]],
                              constant int &nb[[ buffer(7) ]],
                              device const int *b[[ buffer(8) ]],
                              device const int *bo[[ buffer(9) ]],
                              uint tid [[ thread_position_in_grid ]],
                              uint stride [[ threads_per_grid ]]) {

  for (int i = tid; i < nb; i += stride) {
    const int blk_len = b[i];
    const int k = bo[i];
    float tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp += u[gd[k + j] - 1];
    }
    v[dg[k] - 1] = tmp;
  }

  if (o < 0) {
    for (int i = ((-o - 1) + int(tid)); i < m; i += int(stride)) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  } else {
    if ((tid % 2 == 0)) {
      for (int i = ((o - 1) + int(tid)); i < m; i += int(stride)) {
        float tmp = u[gd[i] - 1] + u[gd[i + 1] - 1];
        v[dg[i] - 1] = tmp;
      }
    }
  }
}

/**
 * Gather kernel — multiplication
 */
kernel void gather_kernel_mul(device float *v[[ buffer(0) ]],
                              constant int &m[[ buffer(1) ]],
                              constant int &o[[ buffer(2) ]],
                              device const int *dg[[ buffer(3) ]],
                              device const float *u[[ buffer(4) ]],
                              constant int &n[[ buffer(5) ]],
                              device const int *gd[[ buffer(6) ]],
                              constant int &nb[[ buffer(7) ]],
                              device const int *b[[ buffer(8) ]],
                              device const int *bo[[ buffer(9) ]],
                              uint tid [[ thread_position_in_grid ]],
                              uint stride [[ threads_per_grid ]]) {

  for (int i = tid; i < nb; i += stride) {
    const int blk_len = b[i];
    const int k = bo[i];
    float tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp *= u[gd[k + j] - 1];
    }
    v[dg[k] - 1] = tmp;
  }

  if (o < 0) {
    for (int i = ((-o - 1) + int(tid)); i < m; i += int(stride)) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  } else {
    if ((tid % 2 == 0)) {
      for (int i = ((o - 1) + int(tid)); i < m; i += int(stride)) {
        float tmp = u[gd[i] - 1] * u[gd[i + 1] - 1];
        v[dg[i] - 1] = tmp;
      }
    }
  }
}

/**
 * Gather kernel — minimum
 */
kernel void gather_kernel_min(device float *v[[ buffer(0) ]],
                              constant int &m[[ buffer(1) ]],
                              constant int &o[[ buffer(2) ]],
                              device const int *dg[[ buffer(3) ]],
                              device const float *u[[ buffer(4) ]],
                              constant int &n[[ buffer(5) ]],
                              device const int *gd[[ buffer(6) ]],
                              constant int &nb[[ buffer(7) ]],
                              device const int *b[[ buffer(8) ]],
                              device const int *bo[[ buffer(9) ]],
                              uint tid [[ thread_position_in_grid ]],
                              uint stride [[ threads_per_grid ]]) {

  for (int i = tid; i < nb; i += stride) {
    const int blk_len = b[i];
    const int k = bo[i];
    float tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp = min(u[gd[k + j] - 1], tmp);
    }
    v[dg[k] - 1] = tmp;
  }

  if (o < 0) {
    for (int i = ((-o - 1) + int(tid)); i < m; i += int(stride)) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  } else {
    if ((tid % 2 == 0)) {
      for (int i = ((o - 1) + int(tid)); i < m; i += int(stride)) {
        float tmp = min(u[gd[i] - 1], u[gd[i + 1] - 1]);
        v[dg[i] - 1] = tmp;
      }
    }
  }
}

/**
 * Gather kernel — maximum
 */
kernel void gather_kernel_max(device float *v[[ buffer(0) ]],
                              constant int &m[[ buffer(1) ]],
                              constant int &o[[ buffer(2) ]],
                              device const int *dg[[ buffer(3) ]],
                              device const float *u[[ buffer(4) ]],
                              constant int &n[[ buffer(5) ]],
                              device const int *gd[[ buffer(6) ]],
                              constant int &nb[[ buffer(7) ]],
                              device const int *b[[ buffer(8) ]],
                              device const int *bo[[ buffer(9) ]],
                              uint tid [[ thread_position_in_grid ]],
                              uint stride [[ threads_per_grid ]]) {

  for (int i = tid; i < nb; i += stride) {
    const int blk_len = b[i];
    const int k = bo[i];
    float tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp = max(u[gd[k + j] - 1], tmp);
    }
    v[dg[k] - 1] = tmp;
  }

  if (o < 0) {
    for (int i = ((-o - 1) + int(tid)); i < m; i += int(stride)) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  } else {
    if ((tid % 2 == 0)) {
      for (int i = ((o - 1) + int(tid)); i < m; i += int(stride)) {
        float tmp = max(u[gd[i] - 1], u[gd[i + 1] - 1]);
        v[dg[i] - 1] = tmp;
      }
    }
  }
}

/*
 * Scatter kernel
 */

/**
 * Scatter kernel — copy gathered values back to DOF space
 * u(gd(i)) = v(dg(i))
 */
kernel void scatter_kernel(device const float *v[[ buffer(0) ]],
                           constant int &m[[ buffer(1) ]],
                           device const int *dg[[ buffer(2) ]],
                           device float *u[[ buffer(3) ]],
                           constant int &n[[ buffer(4) ]],
                           device const int *gd[[ buffer(5) ]],
                           constant int &nb[[ buffer(6) ]],
                           device const int *b[[ buffer(7) ]],
                           device const int *bo[[ buffer(8) ]],
                           uint tid [[ thread_position_in_grid ]],
                           uint stride [[ threads_per_grid ]]) {

  for (int i = tid; i < nb; i += stride) {
    const int blk_len = b[i];
    const int k = bo[i];
    float tmp = v[dg[k] - 1];
    for (int j = 0; j < blk_len; j++) {
      u[gd[k + j] - 1] = tmp;
    }
  }

  int facet_offset = 0;
  if (nb > 0) {
    facet_offset = bo[nb - 1] + b[nb - 1];
  }

  for (int i = (facet_offset + int(tid)); i < m; i += int(stride)) {
    u[gd[i] - 1] = v[dg[i] - 1];
  }
}

/*
 * Pack / unpack kernels for MPI communication
 */

/**
 * Pack send buffer from field
 * buf[j] = u[dof[j] - 1]
 */
kernel void gs_pack_kernel(device const float *u[[ buffer(0) ]],
                           device float *buf[[ buffer(1) ]],
                           device const int *dof[[ buffer(2) ]],
                           constant int &n[[ buffer(3) ]],
                           uint tid [[ thread_position_in_grid ]]) {

  if (int(tid) >= n) return;
  buf[tid] = u[dof[tid] - 1];
}

/**
 * Unpack receive buffer — addition
 * Negative dof index means the value must be accumulated atomically.
 */
kernel void gs_unpack_add_kernel(device float *u[[ buffer(0) ]],
                                 device const float *buf[[ buffer(1) ]],
                                 device const int *dof[[ buffer(2) ]],
                                 constant int &n[[ buffer(3) ]],
                                 uint tid [[ thread_position_in_grid ]]) {

  if (int(tid) >= n) return;

  const int idx = dof[tid];
  const float val = buf[tid];
  if (idx < 0) {
    /*
     * Metal supports atomic_fetch_add_explicit on device float via
     * metal::atomic<float> reinterpret. Use atomic_float.
     */
    device atomic_float *addr =
        (device atomic_float *)&u[-idx - 1];
    atomic_fetch_add_explicit(addr, val, memory_order_relaxed);
  } else {
    u[idx - 1] += val;
  }
}

/**
 * Unpack receive buffer — minimum
 * Negative dof index means the value must be reduced atomically.
 */
kernel void gs_unpack_min_kernel(device float *u[[ buffer(0) ]],
                                 device const float *buf[[ buffer(1) ]],
                                 device const int *dof[[ buffer(2) ]],
                                 constant int &n[[ buffer(3) ]],
                                 uint tid [[ thread_position_in_grid ]]) {

  if (int(tid) >= n) return;

  const int idx = dof[tid];
  const float val = buf[tid];
  if (idx < 0) {
    device atomic_float *addr =
        (device atomic_float *)&u[-idx - 1];
    atomic_float_min(addr, val);
  } else {
    u[idx - 1] = min(u[idx - 1], val);
  }
}

/**
 * Unpack receive buffer — maximum
 * Negative dof index means the value must be reduced atomically.
 */
kernel void gs_unpack_max_kernel(device float *u[[ buffer(0) ]],
                                 device const float *buf[[ buffer(1) ]],
                                 device const int *dof[[ buffer(2) ]],
                                 constant int &n[[ buffer(3) ]],
                                 uint tid [[ thread_position_in_grid ]]) {

  if (int(tid) >= n) return;

  const int idx = dof[tid];
  const float val = buf[tid];
  if (idx < 0) {
    device atomic_float *addr =
        (device atomic_float *)&u[-idx - 1];
    atomic_float_max(addr, val);
  } else {
    u[idx - 1] = max(u[idx - 1], val);
  }
}
