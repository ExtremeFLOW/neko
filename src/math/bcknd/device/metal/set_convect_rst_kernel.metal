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
 * Metal compute kernels for the rst convective coefficients (set_convect_rst).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * set_convect_rst — flat per-point, one element per threadgroup
 */

template <int LX>
void set_convect_rst_impl(device float *cr,
                          device float *cs,
                          device float *ct,
                          device const float *cx,
                          device const float *cy,
                          device const float *cz,
                          device const float *drdx,
                          device const float *dsdx,
                          device const float *dtdx,
                          device const float *drdy,
                          device const float *dsdy,
                          device const float *dtdy,
                          device const float *drdz,
                          device const float *dsdz,
                          device const float *dtdz,
                          device const float *w3,
                          uint e, uint iii, uint tpg) {

  const int ele = e * LX * LX * LX;

  for (int ijk = int(iii); ijk < LX * LX * LX; ijk += int(tpg)) {
    const int idx = ijk + ele;
    const float ccx = cx[idx];
    const float ccy = cy[idx];
    const float ccz = cz[idx];
    const float W3 = w3[ijk];

    cr[idx] = W3 * (drdx[idx] * ccx + drdy[idx] * ccy + drdz[idx] * ccz);
    cs[idx] = W3 * (dsdx[idx] * ccx + dsdy[idx] * ccy + dsdz[idx] * ccz);
    ct[idx] = W3 * (dtdx[idx] * ccx + dtdy[idx] * ccy + dtdz[idx] * ccz);
  }
}

#define INSTANTIATE_SET_CONVECT_RST(LX)                                         \
kernel void set_convect_rst_kernel_lx##LX(                                     \
    device float *cr[[ buffer(0) ]],                                           \
    device float *cs[[ buffer(1) ]],                                           \
    device float *ct[[ buffer(2) ]],                                           \
    device const float *cx[[ buffer(3) ]],                                     \
    device const float *cy[[ buffer(4) ]],                                     \
    device const float *cz[[ buffer(5) ]],                                     \
    device const float *drdx[[ buffer(6) ]],                                   \
    device const float *dsdx[[ buffer(7) ]],                                   \
    device const float *dtdx[[ buffer(8) ]],                                   \
    device const float *drdy[[ buffer(9) ]],                                   \
    device const float *dsdy[[ buffer(10) ]],                                  \
    device const float *dtdy[[ buffer(11) ]],                                  \
    device const float *drdz[[ buffer(12) ]],                                  \
    device const float *dsdz[[ buffer(13) ]],                                  \
    device const float *dtdz[[ buffer(14) ]],                                  \
    device const float *w3[[ buffer(15) ]],                                    \
    uint eid [[ threadgroup_position_in_grid ]],                               \
    uint tid [[ thread_index_in_threadgroup ]],                                \
    uint tpg [[ threads_per_threadgroup ]]) {                                  \
  set_convect_rst_impl<LX>(cr, cs, ct, cx, cy, cz,                            \
                           drdx, dsdx, dtdx, drdy, dsdy, dtdy,                 \
                           drdz, dsdz, dtdz, w3, eid, tid, tpg);              \
}

INSTANTIATE_SET_CONVECT_RST(2)
INSTANTIATE_SET_CONVECT_RST(3)
INSTANTIATE_SET_CONVECT_RST(4)
INSTANTIATE_SET_CONVECT_RST(5)
INSTANTIATE_SET_CONVECT_RST(6)
INSTANTIATE_SET_CONVECT_RST(7)
INSTANTIATE_SET_CONVECT_RST(8)
INSTANTIATE_SET_CONVECT_RST(9)
INSTANTIATE_SET_CONVECT_RST(10)
INSTANTIATE_SET_CONVECT_RST(11)
INSTANTIATE_SET_CONVECT_RST(12)
INSTANTIATE_SET_CONVECT_RST(13)
INSTANTIATE_SET_CONVECT_RST(14)
INSTANTIATE_SET_CONVECT_RST(15)
INSTANTIATE_SET_CONVECT_RST(16)
