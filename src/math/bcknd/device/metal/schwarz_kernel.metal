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
 * Metal compute kernels for the additive Schwarz method.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * schwarz_extrude — sum the element "shells", one element per group
 */

template <int NX>
void schwarz_extrude_impl(device float *a1,
                          const int l1,
                          const float f1,
                          device float *a2,
                          const int l2,
                          const float f2,
                          uint eid, uint tid) {

  /*
   * NX == 2 has no interior shell; nothing to extrude. The guarded
   * divisor also avoids a compile-time division-by-zero in that
   * (never-dispatched) template instantiation.
   */
  if (NX <= 2) return;
  const int span = NX - 2;

  const int idx = int(tid);
  const int el = int(eid) * NX * NX * NX;
  const int x = idx % span + 1;
  const int y = idx / span + 1;
  int idx1, idx2;

  idx1 = l1 + x * NX + y * NX * NX + el;
  idx2 = l2 + x * NX + y * NX * NX + el;
  a1[idx1] = f1 * a1[idx1] + f2 * a2[idx2];

  idx1 = (NX - 1 - l1) + x * NX + y * NX * NX + el;
  idx2 = (NX - 1 - l2) + x * NX + y * NX * NX + el;
  a1[idx1] = f1 * a1[idx1] + f2 * a2[idx2];

  threadgroup_barrier(mem_flags::mem_device);

  idx1 = x + l1 * NX + y * NX * NX + el;
  idx2 = x + l2 * NX + y * NX * NX + el;
  a1[idx1] = f1 * a1[idx1] + f2 * a2[idx2];

  idx1 = x + (NX - 1 - l1) * NX + y * NX * NX + el;
  idx2 = x + (NX - 1 - l2) * NX + y * NX * NX + el;
  a1[idx1] = f1 * a1[idx1] + f2 * a2[idx2];

  threadgroup_barrier(mem_flags::mem_device);

  idx1 = x + y * NX + l1 * NX * NX + el;
  idx2 = x + y * NX + l2 * NX * NX + el;
  a1[idx1] = f1 * a1[idx1] + f2 * a2[idx2];

  idx1 = x + y * NX + (NX - 1 - l1) * NX * NX + el;
  idx2 = x + y * NX + (NX - 1 - l2) * NX * NX + el;
  a1[idx1] = f1 * a1[idx1] + f2 * a2[idx2];
}

#define INSTANTIATE_SCHWARZ_EXTRUDE(NX)                                         \
kernel void schwarz_extrude_kernel_nx##NX(                                     \
    device float *a1[[ buffer(0) ]],                                           \
    constant int &l1[[ buffer(1) ]],                                           \
    constant float &f1[[ buffer(2) ]],                                         \
    device float *a2[[ buffer(3) ]],                                           \
    constant int &l2[[ buffer(4) ]],                                           \
    constant float &f2[[ buffer(5) ]],                                         \
    uint eid [[ threadgroup_position_in_grid ]],                               \
    uint tid [[ thread_index_in_threadgroup ]]) {                              \
  schwarz_extrude_impl<NX>(a1, l1, f1, a2, l2, f2, eid, tid);                 \
}

INSTANTIATE_SCHWARZ_EXTRUDE(2)
INSTANTIATE_SCHWARZ_EXTRUDE(3)
INSTANTIATE_SCHWARZ_EXTRUDE(4)
INSTANTIATE_SCHWARZ_EXTRUDE(5)
INSTANTIATE_SCHWARZ_EXTRUDE(6)
INSTANTIATE_SCHWARZ_EXTRUDE(7)
INSTANTIATE_SCHWARZ_EXTRUDE(8)
INSTANTIATE_SCHWARZ_EXTRUDE(9)
INSTANTIATE_SCHWARZ_EXTRUDE(10)
INSTANTIATE_SCHWARZ_EXTRUDE(11)
INSTANTIATE_SCHWARZ_EXTRUDE(12)
INSTANTIATE_SCHWARZ_EXTRUDE(13)
INSTANTIATE_SCHWARZ_EXTRUDE(14)
INSTANTIATE_SCHWARZ_EXTRUDE(15)
INSTANTIATE_SCHWARZ_EXTRUDE(16)

/*
 * schwarz_toext3d — copy regular (nx^3) into extended (nx+2)^3
 */

kernel void schwarz_toext3d_kernel(device float *a[[ buffer(0) ]],
                                   device const float *b[[ buffer(1) ]],
                                   constant int &nx[[ buffer(2) ]],
                                   uint eid [[ threadgroup_position_in_grid ]],
                                   uint tid [[ thread_index_in_threadgroup ]],
                                   uint tpg [[ threads_per_threadgroup ]]) {
  const int idx = int(tid);
  const int str = int(tpg);
  const int nx2 = nx + 2;
  const int el2 = int(eid) * nx2 * nx2 * nx2;
  const int el = int(eid) * nx * nx * nx;

  for (int i = idx; i < nx2 * nx2 * nx2; i += str) {
    a[i + el2] = 0.0f;
  }

  threadgroup_barrier(mem_flags::mem_device);

  for (int ijk = idx; ijk < nx * nx * nx; ijk += str) {
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    a[(i + 1) + (j + 1) * nx2 + (k + 1) * nx2 * nx2 + el2] = b[ijk + el];
  }
}

/*
 * schwarz_toreg3d — copy extended (nx+2)^3 back to regular (nx^3)
 */

kernel void schwarz_toreg3d_kernel(device float *b[[ buffer(0) ]],
                                   device const float *a[[ buffer(1) ]],
                                   constant int &nx[[ buffer(2) ]],
                                   uint eid [[ threadgroup_position_in_grid ]],
                                   uint tid [[ thread_index_in_threadgroup ]],
                                   uint tpg [[ threads_per_threadgroup ]]) {
  const int idx = int(tid);
  const int str = int(tpg);
  const int nx2 = nx + 2;
  const int el2 = int(eid) * nx2 * nx2 * nx2;
  const int el = int(eid) * nx * nx * nx;

  for (int ijk = idx; ijk < nx * nx * nx; ijk += str) {
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    b[ijk + el] = a[(i + 1) + (j + 1) * nx2 + (k + 1) * nx2 * nx2 + el2];
  }
}
