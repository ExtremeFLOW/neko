/*
 Copyright (c) 2022-2025, The Neko Authors
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
 * Device kernel for finding rst via legendre
 */ 

#define DEFINE_FIND_RST_LEGENDRE_KERNEL(LX, CHUNKS)                            \
__kernel                                                                       \
void find_rst_legendre_kernel_lx##LX(__global real * __restrict__ rst,         \
                                     __global const real * pt_x,               \
                                     __global const real * pt_y,               \
                                     __global const real * pt_z,               \
                                     __global const real * x_hat,              \
                                     __global const real * y_hat,              \
                                     __global const real * z_hat,              \
                                     __global real * __restrict__ resx,        \
                                     __global real * __restrict__ resy,        \
                                     __global real * __restrict__ resz,        \
                                     __global const int * el_ids,              \
                                     const int n_pt,                           \
                                     real tol,                                 \
                                     __global real * __restrict__ conv_pts){   \
                                                                               \
  const int pt = get_group_id(0);                                              \
  if (conv_pts[pt] < 0.5) return;                                              \
  const int e = el_ids[pt];                                                    \
  const int elx3 = e*LX*LX*LX;                                                 \
  const int str = get_local_size(1);                                           \
  const int idx = get_local_id(1);                                             \
                                                                               \
  const real one = 1.0;                                                        \
  const real two = 2.0;                                                        \
  __local real dxdr, dydr, dzdr;                                               \
  __local real dxds, dyds, dzds;                                               \
  __local real newx, newy, newz;                                               \
  __local real r_leg[LX];                                                      \
  __local real s_leg[LX];                                                      \
  __local real t_leg[LX];                                                      \
  __local real dr_leg[LX];                                                     \
  __local real ds_leg[LX];                                                     \
  __local real dt_leg[LX];                                                     \
                                                                               \
  __local real xwork[LX*LX];                                                   \
  __local real ywork[LX*LX];                                                   \
  __local real zwork[LX*LX];                                                   \
  __local real xwork2[LX];                                                     \
  __local real ywork2[LX];                                                     \
  __local real zwork2[LX];                                                     \
                                                                               \
  r_leg[0] = 1.0;                                                              \
  s_leg[0] = 1.0;                                                              \
  t_leg[0] = 1.0;                                                              \
  r_leg[1] = rst[3*pt];                                                        \
  s_leg[1] = rst[1+3*pt];                                                      \
  t_leg[1] = rst[2+3*pt];                                                      \
  dr_leg[0] = 0.0;                                                             \
  ds_leg[0] = 0.0;                                                             \
  dt_leg[0] = 0.0;                                                             \
  ds_leg[1] = 1.0;                                                             \
  dr_leg[1] = 1.0;                                                             \
  dt_leg[1] = 1.0;                                                             \
                                                                               \
  for (int ii = 1; ii<LX-1; ii += 1) {                                         \
    real ir = ii;                                                              \
    r_leg[ii+1] =                                                              \
      ((two*ir+one) * rst[3*pt] * r_leg[ii] - ir * r_leg[ii-1] ) / (ir+one);   \
    s_leg[ii+1] =                                                              \
      ((two*ir+one) * rst[3*pt+1] * s_leg[ii] - ir * s_leg[ii-1] ) / (ir+one); \
    t_leg[ii+1] =                                                              \
      ((two*ir+one) * rst[3*pt+2] * t_leg[ii] - ir * t_leg[ii-1] ) / (ir+one); \
    dr_leg[ii+1] = (ir+one) * r_leg[ii] + rst[3*pt]*dr_leg[ii];                \
    ds_leg[ii+1] = (ir+one) * s_leg[ii] + rst[3*pt+1]*ds_leg[ii];              \
    dt_leg[ii+1] = (ir+one) * t_leg[ii] + rst[3*pt+2]*dt_leg[ii];              \
  }                                                                            \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
  for (int ii = idx; ii< LX*LX; ii += str) {                                   \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    int j = ii;                                                                \
    int i = ii - j;                                                            \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];                                   \
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];                                   \
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];                                   \
    }                                                                          \
    xwork[ii] = xtmp;                                                          \
    ywork[ii] = ytmp;                                                          \
    zwork[ii] = ztmp;                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ijk = idx; ijk< LX; ijk += str) {                                   \
    const int jk = ijk;                                                        \
    const int i = ijk - jk;                                                    \
    const int k = jk;                                                          \
    const int j = jk - k;                                                      \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    const int ik2 = i + k*LX;                                                  \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];                                      \
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];                                      \
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];                                      \
    }                                                                          \
    xwork2[ijk] = xtmp;                                                        \
    ywork2[ijk] = ytmp;                                                        \
    zwork2[ijk] = ztmp;                                                        \
  }                                                                            \
                                                                               \
  if(idx==0){                                                                  \
    const int ijk = idx;                                                       \
    const int jk = ijk;                                                        \
    const int i = ijk - jk;                                                    \
    const int k = jk;                                                          \
    const int j = jk - k;                                                      \
    const int ij2 = i + j;                                                     \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];                                   \
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];                                   \
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];                                   \
    }                                                                          \
    newx = xtmp;                                                               \
    newy = ytmp;                                                               \
    newz = ztmp;                                                               \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ii = idx; ii< LX*LX; ii += str) {                                   \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    int j = ii;                                                                \
    int i = ii - j;                                                            \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += dr_leg[i+l]*x_hat[l+LX*j+elx3];                                  \
      ytmp += dr_leg[i+l]*y_hat[l+LX*j+elx3];                                  \
      ztmp += dr_leg[i+l]*z_hat[l+LX*j+elx3];                                  \
    }                                                                          \
    xwork[ii] = xtmp;                                                          \
    ywork[ii] = ytmp;                                                          \
    zwork[ii] = ztmp;                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ijk = idx; ijk< LX; ijk += str) {                                   \
    const int jk = ijk;                                                        \
    const int i = ijk - jk;                                                    \
    const int k = jk;                                                          \
    const int j = jk - k;                                                      \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    const int ik2 = i + k*LX;                                                  \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];                                      \
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];                                      \
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];                                      \
    }                                                                          \
    xwork2[ijk] = xtmp;                                                        \
    ywork2[ijk] = ytmp;                                                        \
    zwork2[ijk] = ztmp;                                                        \
  }                                                                            \
                                                                               \
  if (idx == 0) {                                                              \
    const int ijk = idx;                                                       \
    const int jk = ijk;                                                        \
    const int i = ijk - jk;                                                    \
    const int k = jk ;                                                         \
    const int j = jk - k;                                                      \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    const int ij2 = i + j;                                                     \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];                                   \
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];                                   \
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];                                   \
    }                                                                          \
    dxdr = xtmp;                                                               \
    dydr = ytmp;                                                               \
    dzdr = ztmp;                                                               \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
  for (int ii = idx; ii< LX*LX; ii += str) {                                   \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    int j = ii;                                                                \
    int i = ii - j;                                                            \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];                                   \
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];                                   \
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];                                   \
    }                                                                          \
    xwork[ii] = xtmp;                                                          \
    ywork[ii] = ytmp;                                                          \
    zwork[ii] = ztmp;                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ijk = idx; ijk< LX; ijk += str) {                                   \
    const int jk = ijk;                                                        \
    const int i = ijk - jk;                                                    \
    const int k = jk;                                                          \
    const int j = jk - k;                                                      \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    const int ik2 = i + k*LX;                                                  \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += ds_leg[l+j*LX]*xwork[l+ik2];                                     \
      ytmp += ds_leg[l+j*LX]*ywork[l+ik2];                                     \
      ztmp += ds_leg[l+j*LX]*zwork[l+ik2];                                     \
    }                                                                          \
    xwork2[ijk] = xtmp;                                                        \
    ywork2[ijk] = ytmp;                                                        \
    zwork2[ijk] = ztmp;                                                        \
  }                                                                            \
                                                                               \
  if (idx == 0) {                                                              \
    const int ijk = idx;                                                       \
    const int jk = ijk;                                                        \
    const int i = ijk - jk;                                                    \
    const int k = jk;                                                          \
    const int j = jk - k;                                                      \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    const int ij2 = i + j;                                                     \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];                                   \
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];                                   \
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];                                   \
    }                                                                          \
    dxds = xtmp;                                                               \
    dyds = ytmp;                                                               \
    dzds = ztmp;                                                               \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
  for (int ii = idx; ii< LX*LX; ii += str) {                                   \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    int j = ii;                                                                \
    int i = ii - j;                                                            \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];                                   \
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];                                   \
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];                                   \
    }                                                                          \
    xwork[ii] = xtmp;                                                          \
    ywork[ii] = ytmp;                                                          \
    zwork[ii] = ztmp;                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ijk = idx; ijk< LX; ijk += str) {                                   \
    const int jk = ijk;                                                        \
    const int i = ijk - jk;                                                    \
    const int k = jk;                                                          \
    const int j = jk - k;                                                      \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    const int ik2 = i + k*LX;                                                  \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];                                      \
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];                                      \
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];                                      \
    }                                                                          \
    xwork2[ijk] = xtmp;                                                        \
    ywork2[ijk] = ytmp;                                                        \
    zwork2[ijk] = ztmp;                                                        \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  if( idx == 0){                                                               \
    real xtmp = 0.0;                                                           \
    real ytmp = 0.0;                                                           \
    real ztmp = 0.0;                                                           \
    for( int l = 0; l < LX; l++){                                              \
      xtmp += dt_leg[l]*xwork2[l];                                             \
      ytmp += dt_leg[l]*ywork2[l];                                             \
      ztmp += dt_leg[l]*zwork2[l];                                             \
    }                                                                          \
    real jacdet;                                                               \
    real jacdetinv;                                                            \
    real drdx, drdy, drdz;                                                     \
    real dsdx, dsdy, dsdz;                                                     \
    real dtdx, dtdy, dtdz;                                                     \
    real dxdt, dydt, dzdt;                                                     \
    real rstd[3];                                                              \
    real rstdiff;                                                              \
    real tol2;                                                                 \
    real leg_outside;                                                          \
                                                                               \
    dxdt = xtmp;                                                               \
    dydt = ytmp;                                                               \
    dzdt = ztmp;                                                               \
                                                                               \
    jacdet = (dxdr * dyds * dzdt)                                              \
           + (dxdt * dydr * dzds)                                              \
           + (dxds * dydt * dzdr)                                              \
           - (dxdr * dydt * dzds)                                              \
           - (dxds * dydr * dzdt)                                              \
           - (dxdt * dyds * dzdr);                                             \
                                                                               \
    jacdetinv = one / jacdet;                                                  \
                                                                               \
    drdx =(dyds*dzdt - dydt*dzds);                                             \
    drdy =(dxdt*dzds - dxds*dzdt);                                             \
    drdz =(dxds*dydt - dxdt*dyds);                                             \
    dsdx =(dydt*dzdr - dydr*dzdt);                                             \
    dsdy =(dxdr*dzdt - dxdt*dzdr);                                             \
    dsdz =(dxdt*dydr - dxdr*dydt);                                             \
    dtdx =(dydr*dzds - dyds*dzdr);                                             \
    dtdy =(dxds*dzdr - dxdr*dzds);                                             \
    dtdz =(dxdr*dyds - dxds*dydr);                                             \
                                                                               \
    resx[pt] = pt_x[pt]-newx;                                                  \
    resy[pt] = pt_y[pt]-newy;                                                  \
    resz[pt] = pt_z[pt]-newz;                                                  \
                                                                               \
    rstd[0] = jacdetinv*(drdx*resx[pt]+drdy*resy[pt]+drdz*resz[pt]);           \
    rstd[1] = jacdetinv*(dsdx*resx[pt]+dsdy*resy[pt]+dsdz*resz[pt]);           \
    rstd[2] = jacdetinv*(dtdx*resx[pt]+dtdy*resy[pt]+dtdz*resz[pt]);           \
                                                                               \
    rst[3*pt]     += rstd[0];                                                  \
    rst[3*pt + 1] += rstd[1];                                                  \
    rst[3*pt + 2] += rstd[2];                                                  \
    rstdiff = rstd[0]*rstd[0]+rstd[1]*rstd[1]+rstd[2]*rstd[2];                 \
    conv_pts[pt] = 1.0;                                                        \
    tol2 = tol*tol;                                                            \
    if (rstdiff <= tol2) conv_pts[pt] = 0.0;                                   \
    if (rstdiff > 12.0) conv_pts[pt] = 0.0;                                    \
  }                                                                            \
}

DEFINE_FIND_RST_LEGENDRE_KERNEL(2, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(3, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(4, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(5, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(6, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(7, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(8, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(9, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(10, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(11, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(12, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(13, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(14, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(15, 128)
DEFINE_FIND_RST_LEGENDRE_KERNEL(16, 128)
