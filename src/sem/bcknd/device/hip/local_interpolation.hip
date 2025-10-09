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


#include <stdio.h>
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>

/**
 * Device kernel for finding rst coordinates
 */
template< typename T, typename xT, const int LX, const int CHUNKS>
__global__ void find_rst_legendre_kernel(T * __restrict__ rst,
                                         const T * pt_x,
                                         const T * pt_y,
                                         const T * pt_z,
                                         const T * x_hat,
                                         const T * y_hat,
                                         const T * z_hat,
                                         T * __restrict__ resx,
                                         T * __restrict__ resy,
                                         T * __restrict__ resz,
                                         const int * el_ids,
                                         const int n_pt,
                                         T tol,
                                         T * __restrict__ conv_pts){

  const int pt = blockIdx.x;
  if (conv_pts[pt] < 0.5) return;
  const int e = el_ids[pt];
  const int elx3 = e*LX*LX*LX;
  const int str = blockDim.y;
  const int idx = threadIdx.y;

  const xT one = 1.0;
  const xT two = 2.0;
  __shared__ xT dxdr, dydr, dzdr;
  __shared__ xT dxds, dyds, dzds;
  __shared__ xT newx, newy, newz;
  __shared__ xT r_leg[LX];
  __shared__ xT s_leg[LX];
  __shared__ xT t_leg[LX];
  __shared__ xT dr_leg[LX];
  __shared__ xT ds_leg[LX];
  __shared__ xT dt_leg[LX];


  __shared__ xT xwork[LX*LX];
  __shared__ xT ywork[LX*LX];
  __shared__ xT zwork[LX*LX];
  __shared__ xT xwork2[LX];
  __shared__ xT ywork2[LX];
  __shared__ xT zwork2[LX];

  r_leg[0] = 1.0;
  s_leg[0] = 1.0;
  t_leg[0] = 1.0;
  r_leg[1] = rst[3*pt];
  s_leg[1] = rst[1+3*pt];
  t_leg[1] = rst[2+3*pt];
  dr_leg[0] = 0.0;
  ds_leg[0] = 0.0;
  dt_leg[0] = 0.0;
  ds_leg[1] = 1.0;
  dr_leg[1] = 1.0;
  dt_leg[1] = 1.0;

  for (int ii = 1; ii<LX-1; ii += 1) {
    xT ir = ii;
    r_leg[ii+1] = ((two*ir+one) * rst[3*pt] * r_leg[ii] - ir * r_leg[ii-1] ) / (ir+one);
    s_leg[ii+1] = ((two*ir+one) * rst[3*pt+1] * s_leg[ii] - ir * s_leg[ii-1] ) / (ir+one);
    t_leg[ii+1] = ((two*ir+one) * rst[3*pt+2] * t_leg[ii] - ir * t_leg[ii-1] ) / (ir+one);
    dr_leg[ii+1] = (ir+one) * r_leg[ii] + rst[3*pt]*dr_leg[ii];
    ds_leg[ii+1] = (ir+one) * s_leg[ii] + rst[3*pt+1]*ds_leg[ii];
    dt_leg[ii+1] = (ir+one) * t_leg[ii] + rst[3*pt+2]*dt_leg[ii];
  }
  __syncthreads();
  for (int ii = idx; ii< LX*LX; ii += str) {
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    int j = ii;
    int i = ii - j;
    for( int l = 0; l < LX; l++){
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];
    }
    xwork[ii] = xtmp;
    ywork[ii] = ytmp;
    zwork[ii] = ztmp;
  }

  __syncthreads();

  for (int ijk = idx; ijk< LX; ijk += str) {
    const int jk = ijk;
    const int i = ijk - jk;
    const int k = jk;
    const int j = jk - k;
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    const int ik2 = i + k*LX;
    for( int l = 0; l < LX; l++){
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];
    }
    xwork2[ijk] = xtmp;
    ywork2[ijk] = ytmp;
    zwork2[ijk] = ztmp;
  }

  if(idx==0){
    const int ijk = idx;
    const int jk = ijk;
    const int i = ijk - jk;
    const int k = jk ;
    const int j = jk - k;
    const int ij2 = i + j;
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    for( int l = 0; l < LX; l++){
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];
    }
    newx = xtmp;
    newy = ytmp;
    newz = ztmp;
  //printf("gpu leg and xhat pt rst1 %lf, %lf, %lf, %d, %lf\n",r_leg[5], dr_leg[5],x_hat[elx3],pt, rst[3*pt]);
  }
  __syncthreads();
  for (int ii = idx; ii< LX*LX; ii += str) {
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    int j = ii;
    int i = ii - j;
    for( int l = 0; l < LX; l++){
      xtmp += dr_leg[i+l]*x_hat[l+LX*j+elx3];
      ytmp += dr_leg[i+l]*y_hat[l+LX*j+elx3];
      ztmp += dr_leg[i+l]*z_hat[l+LX*j+elx3];
    }
    xwork[ii] = xtmp;
    ywork[ii] = ytmp;
    zwork[ii] = ztmp;
  }

  __syncthreads();

  for (int ijk = idx; ijk< LX; ijk += str) {
    const int jk = ijk;
    const int i = ijk - jk;
    const int k = jk;
    const int j = jk - k;
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    const int ik2 = i + k*LX;
    for( int l = 0; l < LX; l++){
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];
    }
    xwork2[ijk] = xtmp;
    ywork2[ijk] = ytmp;
    zwork2[ijk] = ztmp;
  }


  if (idx == 0) {
    const int ijk = idx;
    const int jk = ijk;
    const int i = ijk - jk;
    const int k = jk ;
    const int j = jk - k;
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    const int ij2 = i + j;
    for( int l = 0; l < LX; l++){
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];
    }
    dxdr = xtmp;
    dydr = ytmp;
    dzdr = ztmp;
  }

  __syncthreads();
  for (int ii = idx; ii< LX*LX; ii += str) {
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    int j = ii;
    int i = ii - j;
    for( int l = 0; l < LX; l++){
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];
    }
    xwork[ii] = xtmp;
    ywork[ii] = ytmp;
    zwork[ii] = ztmp;
  }

  __syncthreads();

  for (int ijk = idx; ijk< LX; ijk += str) {
    const int jk = ijk;
    const int i = ijk - jk;
    const int k = jk;
    const int j = jk - k;
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    const int ik2 = i + k*LX;
    for( int l = 0; l < LX; l++){
      xtmp += ds_leg[l+j*LX]*xwork[l+ik2];
      ytmp += ds_leg[l+j*LX]*ywork[l+ik2];
      ztmp += ds_leg[l+j*LX]*zwork[l+ik2];
    }
    xwork2[ijk] = xtmp;
    ywork2[ijk] = ytmp;
    zwork2[ijk] = ztmp;
  }


  if (idx == 0) {
    const int ijk = idx;
    const int jk = ijk;
    const int i = ijk - jk;
    const int k = jk ;
    const int j = jk - k;
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    const int ij2 = i + j;
    for( int l = 0; l < LX; l++){
      xtmp += t_leg[l+k*LX]*xwork2[ij2 + l];
      ytmp += t_leg[l+k*LX]*ywork2[ij2 + l];
      ztmp += t_leg[l+k*LX]*zwork2[ij2 + l];
    }
    dxds = xtmp;
    dyds = ytmp;
    dzds = ztmp;
  }

  __syncthreads();
  for (int ii = idx; ii< LX*LX; ii += str) {
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    int j = ii;
    int i = ii - j;
    for( int l = 0; l < LX; l++){
      xtmp += r_leg[i+l]*x_hat[l+LX*j+elx3];
      ytmp += r_leg[i+l]*y_hat[l+LX*j+elx3];
      ztmp += r_leg[i+l]*z_hat[l+LX*j+elx3];
    }
    xwork[ii] = xtmp;
    ywork[ii] = ytmp;
    zwork[ii] = ztmp;
  }

  __syncthreads();

  for (int ijk = idx; ijk< LX; ijk += str) {
    const int jk = ijk;
    const int i = ijk - jk;
    const int k = jk;
    const int j = jk - k;
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    const int ik2 = i + k*LX;
    for( int l = 0; l < LX; l++){
      xtmp += s_leg[l+j*LX]*xwork[l+ik2];
      ytmp += s_leg[l+j*LX]*ywork[l+ik2];
      ztmp += s_leg[l+j*LX]*zwork[l+ik2];
    }
    xwork2[ijk] = xtmp;
    ywork2[ijk] = ytmp;
    zwork2[ijk] = ztmp;
  }

  __syncthreads();

  if( idx == 0){
    xT xtmp = 0.0;
    xT ytmp = 0.0;
    xT ztmp = 0.0;
    for( int l = 0; l < LX; l++){
      xtmp += dt_leg[l]*xwork2[l];
      ytmp += dt_leg[l]*ywork2[l];
      ztmp += dt_leg[l]*zwork2[l];
    }
    xT jacdet;
    xT jacdetinv;
    xT drdx, drdy, drdz;
    xT dsdx, dsdy, dsdz;
    xT dtdx, dtdy, dtdz;
    xT dxdt, dydt, dzdt;
    xT rstd[3];
    xT rstdiff;
    xT tol2;
    xT leg_outside;

    dxdt = xtmp;
    dydt = ytmp;
    dzdt = ztmp;
   //printf("rlegdrleg %lf %lf %f \n",r_leg[5], dr_leg[5],x_hat[elx3]);

    jacdet = (dxdr * dyds * dzdt)
           + (dxdt * dydr * dzds)
           + (dxds * dydt * dzdr)
           - (dxdr * dydt * dzds)
           - (dxds * dydr * dzdt)
           - (dxdt * dyds * dzdr);

    jacdetinv = one / jacdet;
   //check legendre polynoimial, rst, and xhat

    drdx =(dyds*dzdt - dydt*dzds);
    drdy =(dxdt*dzds - dxds*dzdt);
    drdz =(dxds*dydt - dxdt*dyds);
    dsdx =(dydt*dzdr - dydr*dzdt);
    dsdy =(dxdr*dzdt - dxdt*dzdr);
    dsdz =(dxdt*dydr - dxdr*dydt);
    dtdx =(dydr*dzds - dyds*dzdr);
    dtdy =(dxds*dzdr - dxdr*dzds);
    dtdz =(dxdr*dyds - dxds*dydr);
    //printf("newx gpu %lf \n",newx);
    resx[pt] = pt_x[pt]-newx;
    resy[pt] = pt_y[pt]-newy;
    resz[pt] = pt_z[pt]-newz;

    rstd[0] = jacdetinv*(drdx*resx[pt]+drdy*resy[pt]+drdz*resz[pt]);
    rstd[1] = jacdetinv*(dsdx*resx[pt]+dsdy*resy[pt]+dsdz*resz[pt]);
    rstd[2] = jacdetinv*(dtdx*resx[pt]+dtdy*resy[pt]+dtdz*resz[pt]);

    rst[3*pt]     += rstd[0];
    rst[3*pt + 1] += rstd[1];
    rst[3*pt + 2] += rstd[2];
    rstdiff = rstd[0]*rstd[0]+rstd[1]*rstd[1]+rstd[2]*rstd[2];
    conv_pts[pt] = 1.0;
    tol2 = tol*tol;
    if (rstdiff <= tol2) conv_pts[pt] = 0.0;
    if (rstdiff > 12.0) conv_pts[pt] = 0.0;
    //leg_outside = LX*LX*(LX-one)*(LX-one)/4.0;
    //if(r_leg[LX-1]*r_leg[LX-1] > leg_outside) conv_pts[pt] = true;
    //if(s_leg[LX-1]*s_leg[LX-1] > leg_outside) conv_pts[pt] = true;
    //if(t_leg[LX-1]*t_leg[LX-1] > leg_outside) conv_pts[pt] = true;

  }
}


extern "C" {

  /**
   * Fortran wrapper for generating geometric factors
   */
  void hip_find_rst_legendre(void *rst,
                              void *pt_x, void* pt_y, void* pt_z,
                              void *x_hat, void *y_hat, void *z_hat,
                              void *resx, void *resy, void *resz,
                              int *lx, void *el_ids, int *n_pt, real *tol,
                             void *conv_pts) {

    const dim3 nthrds(1, 128, 1);
    const dim3 nblcks((*n_pt), 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

#define RST_CASE(LX)                                                           \
    case LX:                                                                   \
      hipLaunchKernelGGL(                                                      \
         HIP_KERNEL_NAME(find_rst_legendre_kernel<real,real_xp, LX, 128>),     \
         nblcks, nthrds, 0, stream,                                            \
         (real *) rst,(real *) pt_x, (real *) pt_y, (real *) pt_z,             \
         (real *) x_hat, (real *) y_hat, (real *) z_hat,                       \
         (real *) resx, (real *) resy, (real *) resz, (int *) el_ids,          \
         *n_pt, *tol, (real *) conv_pts);                                      \
      HIP_CHECK(hipGetLastError());                                            \
      break

    switch(*lx) {
      RST_CASE(2);
      RST_CASE(3);
      RST_CASE(4);
      RST_CASE(5);
      RST_CASE(6);
      RST_CASE(7);
      RST_CASE(8);
      RST_CASE(9);
      RST_CASE(10);
      RST_CASE(11);
      RST_CASE(12);
      RST_CASE(13);
      RST_CASE(14);
      RST_CASE(15);
      RST_CASE(16);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }
  }
}
