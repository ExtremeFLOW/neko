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
 * Metal compute kernel for finding rst via Legendre polynomials.
 *
 * One threadgroup per point (group size 1 x 128); one Newton step per
 * launch. Unlike the CUDA/OpenCL versions, threadgroup barriers are
 * placed before every thread-0 read of xwork2 — Metal gives no implicit
 * simdgroup lockstep guarantee.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

#define DEFINE_FIND_RST_LEGENDRE_KERNEL(LX)                                    \
kernel void                                                                    \
find_rst_legendre_kernel_lx##LX(                                               \
    device float *rst[[ buffer(0) ]],                                          \
    const device float *pt_x[[ buffer(1) ]],                                   \
    const device float *pt_y[[ buffer(2) ]],                                   \
    const device float *pt_z[[ buffer(3) ]],                                   \
    const device float *x_hat[[ buffer(4) ]],                                  \
    const device float *y_hat[[ buffer(5) ]],                                  \
    const device float *z_hat[[ buffer(6) ]],                                  \
    device float *resx[[ buffer(7) ]],                                         \
    device float *resy[[ buffer(8) ]],                                         \
    device float *resz[[ buffer(9) ]],                                         \
    const device int *el_ids[[ buffer(10) ]],                                  \
    constant int &n_pt[[ buffer(11) ]],                                        \
    constant float &tol[[ buffer(12) ]],                                       \
    device float *conv_pts[[ buffer(13) ]],                                    \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                             \
    uint2 tid [[ thread_position_in_threadgroup ]],                            \
    uint2 tptg [[ threads_per_threadgroup ]])                                  \
{                                                                              \
  const int pt = gid2.x;                                                      \
  if (conv_pts[pt] < 0.5f) return;                                            \
  const int e = el_ids[pt];                                                   \
  const int elx3 = e*LX*LX*LX;                                                \
  const int str = tptg.y;                                                     \
  const int idx = tid.y;                                                      \
                                                                               \
  const float one = 1.0f;                                                     \
  const float two = 2.0f;                                                     \
  threadgroup float dxdr, dydr, dzdr;                                         \
  threadgroup float dxds, dyds, dzds;                                         \
  threadgroup float newx, newy, newz;                                         \
  threadgroup float r_leg[LX];                                                \
  threadgroup float s_leg[LX];                                                \
  threadgroup float t_leg[LX];                                                \
  threadgroup float dr_leg[LX];                                               \
  threadgroup float ds_leg[LX];                                               \
  threadgroup float dt_leg[LX];                                               \
                                                                               \
  threadgroup float xwork[LX*LX];                                             \
  threadgroup float ywork[LX*LX];                                             \
  threadgroup float zwork[LX*LX];                                             \
  threadgroup float xwork2[LX];                                               \
  threadgroup float ywork2[LX];                                               \
  threadgroup float zwork2[LX];                                               \
                                                                               \
  r_leg[0] = 1.0f;                                                            \
  s_leg[0] = 1.0f;                                                            \
  t_leg[0] = 1.0f;                                                            \
  r_leg[1] = rst[3*pt];                                                       \
  s_leg[1] = rst[1+3*pt];                                                     \
  t_leg[1] = rst[2+3*pt];                                                     \
  dr_leg[0] = 0.0f;                                                           \
  ds_leg[0] = 0.0f;                                                           \
  dt_leg[0] = 0.0f;                                                           \
  ds_leg[1] = 1.0f;                                                           \
  dr_leg[1] = 1.0f;                                                           \
  dt_leg[1] = 1.0f;                                                           \
                                                                               \
  for (int ii = 1; ii<LX-1; ii += 1) {                                        \
    float ir = ii;                                                            \
    r_leg[ii+1] =                                                             \
      ((two*ir+one) * rst[3*pt] * r_leg[ii] - ir * r_leg[ii-1] ) / (ir+one);  \
    s_leg[ii+1] =                                                             \
      ((two*ir+one) * rst[3*pt+1] * s_leg[ii] - ir * s_leg[ii-1] ) / (ir+one);\
    t_leg[ii+1] =                                                             \
      ((two*ir+one) * rst[3*pt+2] * t_leg[ii] - ir * t_leg[ii-1] ) / (ir+one);\
    dr_leg[ii+1] = (ir+one) * r_leg[ii] + rst[3*pt]*dr_leg[ii];               \
    ds_leg[ii+1] = (ir+one) * s_leg[ii] + rst[3*pt+1]*ds_leg[ii];             \
    dt_leg[ii+1] = (ir+one) * t_leg[ii] + rst[3*pt+2]*dt_leg[ii];             \
  }                                                                            \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
  for (int ii = idx; ii< LX*LX; ii += str) {                                  \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    int j = ii;                                                               \
    int i = ii - j;                                                           \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];                                  \
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];                                  \
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];                                  \
    }                                                                          \
    xwork[ii] = xtmp;                                                         \
    ywork[ii] = ytmp;                                                         \
    zwork[ii] = ztmp;                                                         \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  for (int ijk = idx; ijk< LX; ijk += str) {                                  \
    const int jk = ijk;                                                       \
    const int i = ijk - jk;                                                   \
    const int k = jk;                                                         \
    const int j = jk - k;                                                     \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    const int ik2 = i + k*LX;                                                 \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];                                     \
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];                                     \
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];                                     \
    }                                                                          \
    xwork2[ijk] = xtmp;                                                       \
    ywork2[ijk] = ytmp;                                                       \
    zwork2[ijk] = ztmp;                                                       \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  if(idx==0){                                                                  \
    const int ijk = idx;                                                      \
    const int jk = ijk;                                                       \
    const int i = ijk - jk;                                                   \
    const int k = jk;                                                         \
    const int j = jk - k;                                                     \
    const int ij2 = i + j;                                                    \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];                                  \
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];                                  \
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];                                  \
    }                                                                          \
    newx = xtmp;                                                              \
    newy = ytmp;                                                              \
    newz = ztmp;                                                              \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  for (int ii = idx; ii< LX*LX; ii += str) {                                  \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    int j = ii;                                                               \
    int i = ii - j;                                                           \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += dr_leg[i+l]*x_hat[l+LX*j+elx3];                                 \
      ytmp += dr_leg[i+l]*y_hat[l+LX*j+elx3];                                 \
      ztmp += dr_leg[i+l]*z_hat[l+LX*j+elx3];                                 \
    }                                                                          \
    xwork[ii] = xtmp;                                                         \
    ywork[ii] = ytmp;                                                         \
    zwork[ii] = ztmp;                                                         \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  for (int ijk = idx; ijk< LX; ijk += str) {                                  \
    const int jk = ijk;                                                       \
    const int i = ijk - jk;                                                   \
    const int k = jk;                                                         \
    const int j = jk - k;                                                     \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    const int ik2 = i + k*LX;                                                 \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];                                     \
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];                                     \
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];                                     \
    }                                                                          \
    xwork2[ijk] = xtmp;                                                       \
    ywork2[ijk] = ytmp;                                                       \
    zwork2[ijk] = ztmp;                                                       \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  if (idx == 0) {                                                             \
    const int ijk = idx;                                                      \
    const int jk = ijk;                                                       \
    const int i = ijk - jk;                                                   \
    const int k = jk ;                                                        \
    const int j = jk - k;                                                     \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    const int ij2 = i + j;                                                    \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];                                  \
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];                                  \
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];                                  \
    }                                                                          \
    dxdr = xtmp;                                                              \
    dydr = ytmp;                                                              \
    dzdr = ztmp;                                                              \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
  for (int ii = idx; ii< LX*LX; ii += str) {                                  \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    int j = ii;                                                               \
    int i = ii - j;                                                           \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];                                  \
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];                                  \
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];                                  \
    }                                                                          \
    xwork[ii] = xtmp;                                                         \
    ywork[ii] = ytmp;                                                         \
    zwork[ii] = ztmp;                                                         \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  for (int ijk = idx; ijk< LX; ijk += str) {                                  \
    const int jk = ijk;                                                       \
    const int i = ijk - jk;                                                   \
    const int k = jk;                                                         \
    const int j = jk - k;                                                     \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    const int ik2 = i + k*LX;                                                 \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += ds_leg[l+j*LX]*xwork[l+ik2];                                    \
      ytmp += ds_leg[l+j*LX]*ywork[l+ik2];                                    \
      ztmp += ds_leg[l+j*LX]*zwork[l+ik2];                                    \
    }                                                                          \
    xwork2[ijk] = xtmp;                                                       \
    ywork2[ijk] = ytmp;                                                       \
    zwork2[ijk] = ztmp;                                                       \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  if (idx == 0) {                                                             \
    const int ijk = idx;                                                      \
    const int jk = ijk;                                                       \
    const int i = ijk - jk;                                                   \
    const int k = jk;                                                         \
    const int j = jk - k;                                                     \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    const int ij2 = i + j;                                                    \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];                                  \
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];                                  \
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];                                  \
    }                                                                          \
    dxds = xtmp;                                                              \
    dyds = ytmp;                                                              \
    dzds = ztmp;                                                              \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
  for (int ii = idx; ii< LX*LX; ii += str) {                                  \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    int j = ii;                                                               \
    int i = ii - j;                                                           \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];                                  \
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];                                  \
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];                                  \
    }                                                                          \
    xwork[ii] = xtmp;                                                         \
    ywork[ii] = ytmp;                                                         \
    zwork[ii] = ztmp;                                                         \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  for (int ijk = idx; ijk< LX; ijk += str) {                                  \
    const int jk = ijk;                                                       \
    const int i = ijk - jk;                                                   \
    const int k = jk;                                                         \
    const int j = jk - k;                                                     \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    const int ik2 = i + k*LX;                                                 \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];                                     \
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];                                     \
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];                                     \
    }                                                                          \
    xwork2[ijk] = xtmp;                                                       \
    ywork2[ijk] = ytmp;                                                       \
    zwork2[ijk] = ztmp;                                                       \
  }                                                                            \
                                                                               \
  threadgroup_barrier(mem_flags::mem_threadgroup);                            \
                                                                               \
  if( idx == 0){                                                              \
    float xtmp = 0.0f;                                                        \
    float ytmp = 0.0f;                                                        \
    float ztmp = 0.0f;                                                        \
    for( int l = 0; l < LX; l++){                                             \
      xtmp += dt_leg[l]*xwork2[l];                                            \
      ytmp += dt_leg[l]*ywork2[l];                                            \
      ztmp += dt_leg[l]*zwork2[l];                                            \
    }                                                                          \
    float jacdet;                                                             \
    float jacdetinv;                                                          \
    float drdx, drdy, drdz;                                                   \
    float dsdx, dsdy, dsdz;                                                   \
    float dtdx, dtdy, dtdz;                                                   \
    float dxdt, dydt, dzdt;                                                   \
    float rstd[3];                                                            \
    float rstdiff;                                                            \
    float tol2;                                                               \
                                                                               \
    dxdt = xtmp;                                                              \
    dydt = ytmp;                                                              \
    dzdt = ztmp;                                                              \
                                                                               \
    jacdet = (dxdr * dyds * dzdt)                                             \
           + (dxdt * dydr * dzds)                                             \
           + (dxds * dydt * dzdr)                                             \
           - (dxdr * dydt * dzds)                                             \
           - (dxds * dydr * dzdt)                                             \
           - (dxdt * dyds * dzdr);                                            \
                                                                               \
    jacdetinv = one / jacdet;                                                 \
                                                                               \
    drdx =(dyds*dzdt - dydt*dzds);                                            \
    drdy =(dxdt*dzds - dxds*dzdt);                                            \
    drdz =(dxds*dydt - dxdt*dyds);                                            \
    dsdx =(dydt*dzdr - dydr*dzdt);                                            \
    dsdy =(dxdr*dzdt - dxdt*dzdr);                                            \
    dsdz =(dxdt*dydr - dxdr*dydt);                                            \
    dtdx =(dydr*dzds - dyds*dzdr);                                            \
    dtdy =(dxds*dzdr - dxdr*dzds);                                            \
    dtdz =(dxdr*dyds - dxds*dydr);                                            \
                                                                               \
    resx[pt] = pt_x[pt]-newx;                                                 \
    resy[pt] = pt_y[pt]-newy;                                                 \
    resz[pt] = pt_z[pt]-newz;                                                 \
                                                                               \
    rstd[0] = jacdetinv*(drdx*resx[pt]+drdy*resy[pt]+drdz*resz[pt]);          \
    rstd[1] = jacdetinv*(dsdx*resx[pt]+dsdy*resy[pt]+dsdz*resz[pt]);          \
    rstd[2] = jacdetinv*(dtdx*resx[pt]+dtdy*resy[pt]+dtdz*resz[pt]);          \
                                                                               \
    rst[3*pt]     += rstd[0];                                                 \
    rst[3*pt + 1] += rstd[1];                                                 \
    rst[3*pt + 2] += rstd[2];                                                 \
    rstdiff = rstd[0]*rstd[0]+rstd[1]*rstd[1]+rstd[2]*rstd[2];                \
    conv_pts[pt] = 1.0f;                                                      \
    tol2 = tol*tol;                                                           \
    if (rstdiff <= tol2) conv_pts[pt] = 0.0f;                                 \
    if (rstdiff > 12.0f) conv_pts[pt] = 0.0f;                                 \
  }                                                                            \
}

DEFINE_FIND_RST_LEGENDRE_KERNEL(2)
DEFINE_FIND_RST_LEGENDRE_KERNEL(3)
DEFINE_FIND_RST_LEGENDRE_KERNEL(4)
DEFINE_FIND_RST_LEGENDRE_KERNEL(5)
DEFINE_FIND_RST_LEGENDRE_KERNEL(6)
DEFINE_FIND_RST_LEGENDRE_KERNEL(7)
DEFINE_FIND_RST_LEGENDRE_KERNEL(8)
DEFINE_FIND_RST_LEGENDRE_KERNEL(9)
DEFINE_FIND_RST_LEGENDRE_KERNEL(10)
DEFINE_FIND_RST_LEGENDRE_KERNEL(11)
DEFINE_FIND_RST_LEGENDRE_KERNEL(12)
DEFINE_FIND_RST_LEGENDRE_KERNEL(13)
DEFINE_FIND_RST_LEGENDRE_KERNEL(14)
DEFINE_FIND_RST_LEGENDRE_KERNEL(15)
DEFINE_FIND_RST_LEGENDRE_KERNEL(16)
