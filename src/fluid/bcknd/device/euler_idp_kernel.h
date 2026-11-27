/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of the authors nor the names of its contributors may be
     used to endorse or promote products derived from this software without
     specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
*/

#include <cmath>
#include <cfloat>
#include <device/device_config.h>

namespace {

template<typename T> __device__ __forceinline__ T rmin(T a, T b) {
  return a < b ? a : b;
}

template<typename T> __device__ __forceinline__ T rmax(T a, T b) {
  return a > b ? a : b;
}

template<typename T> __device__ __forceinline__ T real_epsilon();
template<typename T> __device__ __forceinline__ T real_tiny();
template<typename T> __device__ __forceinline__ T real_huge();

template<> __device__ __forceinline__ float real_epsilon<float>() {
  return FLT_EPSILON;
}

template<> __device__ __forceinline__ double real_epsilon<double>() {
  return DBL_EPSILON;
}

template<> __device__ __forceinline__ float real_tiny<float>() {
  return FLT_MIN;
}

template<> __device__ __forceinline__ double real_tiny<double>() {
  return DBL_MIN;
}

template<> __device__ __forceinline__ float real_huge<float>() {
  return FLT_MAX;
}

template<> __device__ __forceinline__ double real_huge<double>() {
  return DBL_MAX;
}

template<typename T> __device__ __forceinline__ T atomic_min_real(T *p, T v);
template<typename T> __device__ __forceinline__ T atomic_max_real(T *p, T v);

template<> __device__ __forceinline__ float atomic_min_real(float *p, float v) {
  int *address = reinterpret_cast<int *>(p);
  int old = *address;
  while (v < __int_as_float(old)) {
    const int assumed = old;
    old = atomicCAS(address, assumed, __float_as_int(v));
    if (old == assumed) break;
  }
  return __int_as_float(old);
}

template<> __device__ __forceinline__ float atomic_max_real(float *p, float v) {
  int *address = reinterpret_cast<int *>(p);
  int old = *address;
  while (v > __int_as_float(old)) {
    const int assumed = old;
    old = atomicCAS(address, assumed, __float_as_int(v));
    if (old == assumed) break;
  }
  return __int_as_float(old);
}

template<> __device__ __forceinline__ double atomic_min_real(double *p,
                                                               double v) {
  unsigned long long *address = reinterpret_cast<unsigned long long *>(p);
  unsigned long long old = *address;
  while (v < __longlong_as_double(old)) {
    const unsigned long long assumed = old;
    old = atomicCAS(address, assumed, __double_as_longlong(v));
    if (old == assumed) break;
  }
  return __longlong_as_double(old);
}

template<> __device__ __forceinline__ double atomic_max_real(double *p,
                                                               double v) {
  unsigned long long *address = reinterpret_cast<unsigned long long *>(p);
  unsigned long long old = *address;
  while (v > __longlong_as_double(old)) {
    const unsigned long long assumed = old;
    old = atomicCAS(address, assumed, __double_as_longlong(v));
    if (old == assumed) break;
  }
  return __longlong_as_double(old);
}

template<typename T>
__device__ __forceinline__ void load_state(const T *rho, const T *mx,
                                            const T *my, const T *mz,
                                            const T *energy, int i, T s[5]) {
  s[0] = rho[i]; s[1] = mx[i]; s[2] = my[i]; s[3] = mz[i]; s[4] = energy[i];
}

template<typename T>
__device__ __forceinline__ T internal_energy(const T s[5]) {
  return s[4] - T(0.5) *
    (s[1]*s[1] + s[2]*s[2] + s[3]*s[3]) / s[0];
}

template<typename T>
__device__ __forceinline__ bool state_is_admissible(const T s[5], T floor) {
  for (int c = 0; c < 5; ++c) if (!isfinite(s[c])) return false;
  if (s[0] <= T(0)) return false;
  const T margin = s[0] * (s[4] - floor) - T(0.5) *
    (s[1]*s[1] + s[2]*s[2] + s[3]*s[3]);
  return isfinite(margin) && margin >= T(0);
}

template<typename T>
__device__ __forceinline__ T specific_entropy(const T s[5], T gamma) {
  const T huge = real_huge<T>();
  for (int c = 0; c < 5; ++c) if (!isfinite(s[c])) return -huge;
  if (s[0] <= T(0)) return -huge;
  const T pressure = (gamma - T(1)) * internal_energy(s);
  if (!isfinite(pressure) || pressure <= T(0)) return -huge;
  return log(pressure) - gamma * log(s[0]);
}

template<typename T>
__device__ __forceinline__ T entropy_tolerance(const T s[5], T lower) {
  const T eps = real_epsilon<T>();
  const T tiny = real_tiny<T>();
  T conditioning = T(1);
  bool finite = s[0] > T(0);
  for (int c = 0; c < 5; ++c) finite = finite && isfinite(s[c]);
  if (finite) {
    const T kinetic = T(0.5) *
      (s[1]*s[1] + s[2]*s[2] + s[3]*s[3]) / s[0];
    const T internal = s[4] - kinetic;
    if (isfinite(internal) && internal > tiny)
      conditioning = (fabs(s[4]) + fabs(kinetic)) / internal;
  }
  return rmin(sqrt(eps), T(64) * eps *
              rmax(T(1), rmax(fabs(lower), conditioning)));
}

template<typename T>
__device__ __forceinline__ bool entropy_is_admissible(const T s[5], T gamma,
                                                       T lower, T floor) {
  if (!state_is_admissible(s, floor)) return false;
  const T entropy = specific_entropy(s, gamma);
  return isfinite(entropy) && entropy >= lower - entropy_tolerance(s, lower);
}

template<typename T>
__device__ __forceinline__ void fluxes(const T s[5], T gamma,
                                        T fx[5], T fy[5], T fz[5]) {
  const T inv_rho = T(1) / s[0];
  const T u = s[1] * inv_rho;
  const T v = s[2] * inv_rho;
  const T w = s[3] * inv_rho;
  const T p = (gamma - T(1)) * internal_energy(s);
  fx[0] = s[1]; fy[0] = s[2]; fz[0] = s[3];
  fx[1] = s[1]*u + p; fy[1] = s[1]*v;     fz[1] = s[1]*w;
  fx[2] = s[2]*u;     fy[2] = s[2]*v + p; fz[2] = s[2]*w;
  fx[3] = s[3]*u;     fy[3] = s[3]*v;     fz[3] = s[3]*w + p;
  fx[4] = (s[4]+p)*u; fy[4] = (s[4]+p)*v; fz[4] = (s[4]+p)*w;
}

template<typename T>
__device__ __forceinline__ T ordered_wave_speed(const T left[5],
                                                 const T right[5],
                                                 const T normal[3], T gamma) {
  const T pl = (gamma - T(1)) * internal_energy(left);
  const T pr = (gamma - T(1)) * internal_energy(right);
  const T al = sqrt(gamma * pl / left[0]);
  const T ar = sqrt(gamma * pr / right[0]);
  const T ul = (left[1]*normal[0] + left[2]*normal[1] +
                left[3]*normal[2]) / left[0];
  const T ur = (right[1]*normal[0] + right[2]*normal[1] +
                right[3]*normal[2]) / right[0];
  T pmin, pmax, dmin, amin, amax;
  if (pl <= pr) { pmin=pl; pmax=pr; dmin=left[0]; amin=al; amax=ar; }
  else { pmin=pr; pmax=pl; dmin=right[0]; amin=ar; amax=al; }
  const T exponent = (gamma-T(1))/(T(2)*gamma);
  const T ratio = pow(pmin/pmax, exponent);
  const T phi_min = T(2)*amax*(ratio-T(1))/(gamma-T(1)) + ur-ul;
  if (phi_min >= T(0))
    return rmax(rmax(-(ul-al), T(0)), rmax(ur+ar, T(0)));
  const T ca = T(2)/((gamma+T(1))*dmin);
  const T cb = pmin*(gamma-T(1))/(gamma+T(1));
  const T phi_max = (pmax-pmin)*sqrt(ca/(pmax+cb)) + ur-ul;
  const T numerator = amin+amax-T(0.5)*(gamma-T(1))*(ur-ul);
  const T p_two = pmin * pow(numerator/(amin+amax*ratio),
                              T(2)*gamma/(gamma-T(1)));
  const T p_upper = phi_max < T(0) ? p_two : rmin(pmax, p_two);
  const T sl = ul - al*sqrt(T(1) + rmax((p_upper-pl)/pl, T(0)) *
                                (gamma+T(1))/(T(2)*gamma));
  const T sr = ur + ar*sqrt(T(1) + rmax((p_upper-pr)/pr, T(0)) *
                                (gamma+T(1))/(T(2)*gamma));
  return rmax(rmax(-sl, T(0)), rmax(sr, T(0)));
}

template<typename T>
__device__ __forceinline__ T maximum_wave_speed(const T left[5],
                                                 const T right[5],
                                                 const T normal[3], T gamma) {
  T reverse[3] = {-normal[0], -normal[1], -normal[2]};
  return rmax(ordered_wave_speed(left, right, normal, gamma),
              ordered_wave_speed(right, left, reverse, gamma));
}

template<typename T>
__device__ __forceinline__ T energy_margin(const T s[5], T floor) {
  return s[4] - floor - T(0.5) *
    (s[1]*s[1] + s[2]*s[2] + s[3]*s[3]) / s[0];
}

template<typename T>
__device__ __forceinline__ T energy_margin_derivative(const T s[5],
                                                       const T q[5]) {
  const T m2 = s[1]*s[1] + s[2]*s[2] + s[3]*s[3];
  return q[4] - (s[1]*q[1]+s[2]*q[2]+s[3]*q[3])/s[0] +
    T(0.5)*m2*q[0]/(s[0]*s[0]);
}

template<typename T>
__device__ __forceinline__ T entropy_margin(const T s[5], T lower, T gamma) {
  return (gamma-T(1))*energy_margin(s, T(0)) - exp(lower)*pow(s[0], gamma);
}

template<typename T>
__device__ __forceinline__ T entropy_margin_derivative(const T s[5],
                                                        const T q[5],
                                                        T lower, T gamma) {
  return (gamma-T(1))*energy_margin_derivative(s,q) -
    exp(lower)*gamma*pow(s[0],gamma-T(1))*q[0];
}

template<typename T>
__device__ void limit_energy(const T base[5], const T correction[5], T floor,
                             T &limit) {
  const T eps = real_epsilon<T>();
  T left=T(0), right=limit, trial[5];
  T vl=energy_margin(base,floor);
  for (int c=0;c<5;++c) trial[c]=base[c]+right*correction[c];
  T vr=energy_margin(trial,floor);
  const T tol=T(256)*eps*rmax(T(1),rmax(fabs(vl),fabs(vr)));
  if (vl <= T(0)) { limit=T(0); return; }
  for (int it=0; it<16; ++it) {
    if (right-left <= T(256)*eps*rmax(T(1),right) || vl <= vr) break;
    const T slope=(vr-vl)/(right-left);
    T next=left-vl/slope;
    if (!isfinite(next) || next<=left || next>=right) break;
    for (int c=0;c<5;++c) trial[c]=base[c]+next*correction[c];
    T vn=energy_margin(trial,floor);
    if (!isfinite(vn) || vn < -tol) break;
    left=next; vl=vn;
    if (fabs(vn)<=tol) break;
    for (int c=0;c<5;++c) trial[c]=base[c]+right*correction[c];
    const T ds=energy_margin_derivative(trial,correction);
    if (!isfinite(ds) || ds>=T(0)) break;
    next=right-vr/ds;
    if (!isfinite(next) || next<=left || next>=right) break;
    for (int c=0;c<5;++c) trial[c]=base[c]+next*correction[c];
    vn=energy_margin(trial,floor);
    if (!isfinite(vn) || vn>tol) break;
    right=next; vr=vn;
    if (fabs(vn)<=tol) break;
  }
  limit=left;
}

template<typename T>
__device__ void limit_entropy(const T base[5], const T correction[5], T lower,
                              T gamma, T &limit) {
  const T eps=real_epsilon<T>();
  T left=T(0), right=limit, trial[5];
  T vl=entropy_margin(base,lower,gamma);
  for (int c=0;c<5;++c) trial[c]=base[c]+right*correction[c];
  T vr=entropy_margin(trial,lower,gamma);
  const T tol=T(256)*eps*rmax(T(1),rmax(fabs(vl),fabs(vr)));
  if (vl<=T(0)) { limit=T(0); return; }
  for (int it=0;it<16;++it) {
    if (right-left<=T(256)*eps*rmax(T(1),right) || vl<=vr) break;
    const T slope=(vr-vl)/(right-left);
    T next=left-vl/slope;
    if (!isfinite(next)||next<=left||next>=right) break;
    for(int c=0;c<5;++c) trial[c]=base[c]+next*correction[c];
    T vn=entropy_margin(trial,lower,gamma);
    if(!isfinite(vn)||vn < -tol) break;
    left=next; vl=vn;
    if(fabs(vn)<=tol) break;
    for(int c=0;c<5;++c) trial[c]=base[c]+right*correction[c];
    const T ds=entropy_margin_derivative(trial,correction,lower,gamma);
    if(!isfinite(ds)||ds>=T(0)) break;
    next=right-vr/ds;
    if(!isfinite(next)||next<=left||next>=right) break;
    for(int c=0;c<5;++c) trial[c]=base[c]+next*correction[c];
    vn=entropy_margin(trial,lower,gamma);
    if(!isfinite(vn)||vn>tol) break;
    right=next; vr=vn;
    if(fabs(vn)<=tol) break;
  }
  limit=left;
}

template<typename T>
__device__ void entropy_bisection(const T base[5], const T correction[5],
                                  T lower_entropy, T gamma, T floor,
                                  T &limit) {
  T trial[5];
  for(int c=0;c<5;++c) trial[c]=base[c]+limit*correction[c];
  if(entropy_is_admissible(trial,gamma,lower_entropy,floor)) return;
  if(!entropy_is_admissible(base,gamma,lower_entropy,floor)) {
    limit=T(0); return;
  }
  T left=T(0), right=limit;
  for(int it=0;it<64;++it) {
    const T mid=left+T(0.5)*(right-left);
    if(mid<=left||mid>=right) break;
    for(int c=0;c<5;++c) trial[c]=base[c]+mid*correction[c];
    if(entropy_is_admissible(trial,gamma,lower_entropy,floor)) left=mid;
    else right=mid;
  }
  limit=left;
}

template<typename T>
__device__ T limit_endpoint(const T base[5], const T correction[5], T dl,
                            T du, T sl, T gamma, T floor, bool enforce_energy,
                            bool enforce_entropy, bool check_base,
                            bool &density_limited, bool &energy_limited,
                            bool &entropy_limited) {
  const T eps=real_epsilon<T>();
  T limit=T(1), trial[5];
  density_limited=false; energy_limited=false; entropy_limited=false;
  if(check_base) {
    bool finite=true;
    for(int c=0;c<5;++c) finite=finite&&isfinite(base[c]);
    if(!finite) { energy_limited=enforce_energy;
      entropy_limited=enforce_entropy; return T(0); }
    if(base[0]<=T(0)||base[0]<dl||base[0]>du) {
      density_limited=true; return T(0); }
    if(enforce_energy&&!state_is_admissible(base,floor)) {
      energy_limited=true; return T(0); }
    if(enforce_entropy&&!entropy_is_admissible(base,gamma,sl,floor)) {
      entropy_limited=true; return T(0); }
  }
  if(correction[0]<T(0)) limit=rmin(limit,(base[0]-dl)/(-correction[0]));
  else if(correction[0]>T(0)) limit=rmin(limit,(du-base[0])/correction[0]);
  if(correction[0]!=T(0)) {
    limit=rmax(T(0),limit); density_limited=limit<T(1);
    if(density_limited&&limit>T(0)) limit*=T(1)-T(32)*eps;
  }
  for(int c=0;c<5;++c) trial[c]=base[c]+limit*correction[c];
  energy_limited=enforce_energy&&!state_is_admissible(trial,floor);
  entropy_limited=enforce_entropy&&
    !entropy_is_admissible(trial,gamma,sl,floor);
  if(!energy_limited&&!entropy_limited) return limit;
  if(energy_limited) limit_energy(base,correction,floor,limit);
  for(int c=0;c<5;++c) trial[c]=base[c]+limit*correction[c];
  if(enforce_entropy&&!entropy_is_admissible(trial,gamma,sl,floor)) {
    limit_entropy(base,correction,sl,gamma,limit); entropy_limited=true;
  }
  if(limit>T(0)) limit*=T(1)-T(32)*eps;
  if(enforce_entropy&&limit>T(0))
    entropy_bisection(base,correction,sl,gamma,floor,limit);
  return limit;
}

template<typename T>
__device__ T floor_timestep(const T state[5], const T residual[5], T floor,
                            T upper) {
  const T eps=real_epsilon<T>();
  T limit=rmax(T(0),upper), trial[5];
  if(limit<=T(0)) return limit;
  if(residual[0]>T(0)) {
    const T density_limit=state[0]/residual[0];
    if(density_limit<=limit)
      limit=rmax(T(0),density_limit*(T(1)-T(64)*eps));
  }
  if(limit<=T(0)) return limit;
  for(int c=0;c<5;++c) trial[c]=state[c]-limit*residual[c];
  if(state_is_admissible(trial,floor)) return limit;
  if(!state_is_admissible(state,floor)) return T(0);
  // Leave room for roundoff when the caller forms the candidate again.
  const T guarded_floor = floor + T(16)*eps*
    rmax(fabs(state[4]), limit*fabs(residual[4]));
  T left=T(0), right=limit;
  for(int it=0;it<64;++it) {
    const T mid=T(0.5)*(left+right);
    for(int c=0;c<5;++c) trial[c]=state[c]-mid*residual[c];
    if(state_is_admissible(trial,guarded_floor)) left=mid; else right=mid;
  }
  return left;
}

template<typename T>
__global__ void primitives_kernel(const T *rho,const T *mx,const T *my,
  const T *mz,const T *energy,T *u,T *v,T *w,T *p,T *sound,T *internal,
  T *status,T gamma,T floor,int n) {
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x) {
    T s[5]; load_state(rho,mx,my,mz,energy,i,s);
    status[i]=T(0); u[i]=v[i]=w[i]=p[i]=sound[i]=internal[i]=T(0);
    if(!isfinite(gamma)||gamma<=T(1)||!state_is_admissible(s,floor)) {
      status[i]=T(1); continue;
    }
    u[i]=s[1]/s[0]; v[i]=s[2]/s[0]; w[i]=s[3]/s[0];
    internal[i]=internal_energy(s); p[i]=(gamma-T(1))*internal[i];
    sound[i]=sqrt(gamma*p[i]/s[0]);
    if(!isfinite(sound[i])) status[i]=T(1);
  }
}

template<typename T>
__global__ void graph_viscosity_kernel(const T *rho,const T *mx,const T *my,
  const T *mz,const T *energy,const T *user_speed,int has_user,
  const int *left,const int *right,const T *coef,T *edge_visc,T *vis_sum,
  T *edge_speed,T gamma,int nedge) {
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nedge;e+=blockDim.x*gridDim.x) {
    const int a=left[e], b=right[e];
    const T cx=coef[3*e],cy=coef[3*e+1],cz=coef[3*e+2];
    const T norm=sqrt(cx*cx+cy*cy+cz*cz);
    T wave;
    if(has_user) wave=rmax(user_speed[a],user_speed[b]);
    else {
      T l[5],r[5],normal[3]={cx/norm,cy/norm,cz/norm};
      load_state(rho,mx,my,mz,energy,a,l); load_state(rho,mx,my,mz,energy,b,r);
      wave=maximum_wave_speed(l,r,normal,gamma);
    }
    edge_speed[e]=wave; edge_visc[e]=norm*wave;
    atomicAdd(vis_sum+a,edge_visc[e]); atomicAdd(vis_sum+b,edge_visc[e]);
  }
}

template<typename T>
__global__ void graph_rate_kernel(const T *vis,const T *mass,T *rate,int n) {
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x)
    rate[i]=T(2)*vis[i]/mass[i];
}

template<typename T>
__global__ void nodal_wave_kernel(const T *u,const T *v,const T *w,
  const T *sound,T *wave,int n) {
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x)
    wave[i]=sqrt(u[i]*u[i]+v[i]*v[i]+w[i]*w[i])+sound[i];
}

template<typename T>
__global__ void low_init_kernel(const T *rho,const T *mx,const T *my,
  const T *mz,const T *energy,const T *diag,T *r0,T *r1,T *r2,T *r3,T *r4,
  T gamma,int periodic,int n) {
  T *out[5]={r0,r1,r2,r3,r4};
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x) {
    if(periodic) { for(int c=0;c<5;++c) out[c][i]=T(0); continue; }
    T s[5],fx[5],fy[5],fz[5]; load_state(rho,mx,my,mz,energy,i,s);
    fluxes(s,gamma,fx,fy,fz);
    for(int c=0;c<5;++c)
      out[c][i]=diag[3*i]*fx[c]+diag[3*i+1]*fy[c]+diag[3*i+2]*fz[c];
  }
}

template<typename T>
__global__ void low_edge_kernel(const T *rho,const T *mx,const T *my,
  const T *mz,const T *energy,const int *left,const int *right,const T *coef,
  const T *visc,T *r0,T *r1,T *r2,T *r3,T *r4,T gamma,int periodic,int nedge) {
  T *out[5]={r0,r1,r2,r3,r4};
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nedge;e+=blockDim.x*gridDim.x) {
    const int a=left[e],b=right[e]; T sa[5],sb[5],fax[5],fay[5],faz[5];
    T fbx[5],fby[5],fbz[5]; load_state(rho,mx,my,mz,energy,a,sa);
    load_state(rho,mx,my,mz,energy,b,sb); fluxes(sa,gamma,fax,fay,faz);
    fluxes(sb,gamma,fbx,fby,fbz);
    const T cx=coef[3*e],cy=coef[3*e+1],cz=coef[3*e+2];
    for(int c=0;c<5;++c) {
      const T fa=cx*fax[c]+cy*fay[c]+cz*faz[c];
      const T fb=cx*fbx[c]+cy*fby[c]+cz*fbz[c];
      const T d=visc[e]*(sb[c]-sa[c]);
      if(periodic) { atomicAdd(out[c]+a,fb-fa-d); atomicAdd(out[c]+b,fb-fa+d); }
      else { atomicAdd(out[c]+a,fb-d); atomicAdd(out[c]+b,-fa+d); }
    }
  }
}

template<typename T>
__global__ void scale_residual_kernel(T *r0,T *r1,T *r2,T *r3,T *r4,
                                      const T *binv,int n) {
  T *out[5]={r0,r1,r2,r3,r4};
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x)
    for(int c=0;c<5;++c) out[c][i]*=binv[i];
}

template<typename T>
__global__ void correction_residual_init_kernel(const T *rho,const T *mx,
  const T *my,const T *mz,const T *energy,const T *diag,
  const T *drdx,const T *drdy,const T *drdz,const T *dsdx,const T *dsdy,
  const T *dsdz,const T *dtdx,const T *dtdy,const T *dtdz,
  const T *dx,const T *dy,const T *dz,const T *wx,const T *wy,const T *wz,
  T *r0,T *r1,T *r2,T *r3,T *r4,T gamma,int lx,int ly,int lz,int n) {
  T *out[5]={r0,r1,r2,r3,r4};
  for(int node=blockIdx.x*blockDim.x+threadIdx.x;node<n;
      node+=blockDim.x*gridDim.x) {
    int q=node; const int i=q%lx; q/=lx; const int j=q%ly; q/=ly;
    const int k=q%lz; const int e=q/lz;
    T s[5],fx[5],fy[5],fz[5]; load_state(rho,mx,my,mz,energy,node,s);
    fluxes(s,gamma,fx,fy,fz); T high[5]={0,0,0,0,0};
    const T weight=wx[i]*wy[j]*wz[k];
    for(int l=0;l<lx;++l) {
      const int nl=((e*lz+k)*ly+j)*lx+l; T sl[5],x[5],y[5],z[5];
      load_state(rho,mx,my,mz,energy,nl,sl); fluxes(sl,gamma,x,y,z);
      const T d=weight*dx[i+l*lx];
      for(int c=0;c<5;++c)
        high[c]+=d*(drdx[nl]*x[c]+drdy[nl]*y[c]+drdz[nl]*z[c]);
    }
    for(int l=0;l<ly;++l) {
      const int nl=((e*lz+k)*ly+l)*lx+i; T sl[5],x[5],y[5],z[5];
      load_state(rho,mx,my,mz,energy,nl,sl); fluxes(sl,gamma,x,y,z);
      const T d=weight*dy[j+l*ly];
      for(int c=0;c<5;++c)
        high[c]+=d*(dsdx[nl]*x[c]+dsdy[nl]*y[c]+dsdz[nl]*z[c]);
    }
    if(lz>1) for(int l=0;l<lz;++l) {
      const int nl=((e*lz+l)*ly+j)*lx+i; T sl[5],x[5],y[5],z[5];
      load_state(rho,mx,my,mz,energy,nl,sl); fluxes(sl,gamma,x,y,z);
      const T d=weight*dz[k+l*lz];
      for(int c=0;c<5;++c)
        high[c]+=d*(dtdx[nl]*x[c]+dtdy[nl]*y[c]+dtdz[nl]*z[c]);
    }
    for(int c=0;c<5;++c)
      out[c][node]=diag[3*node]*fx[c]+diag[3*node+1]*fy[c]+
        diag[3*node+2]*fz[c]-high[c];
  }
}

template<typename T>
__global__ void correction_residual_edge_kernel(const T *rho,const T *mx,
  const T *my,const T *mz,const T *energy,const int *left,const int *right,
  const T *coef,T *r0,T *r1,T *r2,T *r3,T *r4,T gamma,int nedge) {
  T *out[5]={r0,r1,r2,r3,r4};
  for(int edge=blockIdx.x*blockDim.x+threadIdx.x;edge<nedge;
      edge+=blockDim.x*gridDim.x) {
    const int a=left[edge],b=right[edge]; T sa[5],sb[5],fax[5],fay[5],faz[5];
    T fbx[5],fby[5],fbz[5]; load_state(rho,mx,my,mz,energy,a,sa);
    load_state(rho,mx,my,mz,energy,b,sb); fluxes(sa,gamma,fax,fay,faz);
    fluxes(sb,gamma,fbx,fby,fbz); const T cx=coef[3*edge];
    const T cy=coef[3*edge+1],cz=coef[3*edge+2];
    for(int c=0;c<5;++c) {
      atomicAdd(out[c]+a,cx*fbx[c]+cy*fby[c]+cz*fbz[c]);
      atomicAdd(out[c]+b,-cx*fax[c]-cy*fay[c]-cz*faz[c]);
    }
  }
}

template<typename T>
__global__ void reconstruct_general_kernel(const T *r0,const T *r1,
  const T *r2,const T *r3,const T *r4,const T *wx,const T *wy,const T *wz,
  T *correction,T *error,int lx,int ly,int lz,int edges_per_element,int nelv) {
  const T *res[5]={r0,r1,r2,r3,r4};
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nelv;e+=blockDim.x*gridDim.x) {
    T sx=T(0),sy=T(0),sz=T(0); for(int i=0;i<lx;++i)sx+=wx[i];
    for(int j=0;j<ly;++j)sy+=wy[j]; for(int k=0;k<lz;++k)sz+=wz[k];
    T element_sum[5]={0,0,0,0,0},scale[5]={0,0,0,0,0};
    for(int k=0;k<lz;++k)for(int j=0;j<ly;++j)for(int i=0;i<lx;++i) {
      const int node=((e*lz+k)*ly+j)*lx+i;
      for(int c=0;c<5;++c){element_sum[c]+=res[c][node];scale[c]+=fabs(res[c][node]);}
    }
    for(int c=0;c<5;++c)
      atomic_max_real(error+c,fabs(element_sum[c])/rmax(T(1),scale[c]));
    const int base=e*edges_per_element;
    int edge=base;
    for(int k=0;k<lz;++k)for(int j=0;j<ly;++j) {
      T line[5]={0,0,0,0,0},cum[5]={0,0,0,0,0};
      for(int i=0;i<lx;++i){const int n=((e*lz+k)*ly+j)*lx+i;
        for(int c=0;c<5;++c)line[c]+=res[c][n];}
      for(int i=0;i<lx-1;++i){const int n=((e*lz+k)*ly+j)*lx+i;
        for(int c=0;c<5;++c){cum[c]+=res[c][n]-wx[i]/sx*line[c];
          correction[5*edge+c]=cum[c];} ++edge;}
    }
    const int ybase=base+(lx-1)*ly*lz; edge=ybase;
    for(int k=0;k<lz;++k) {
      T plane[5]={0,0,0,0,0};
      for(int j=0;j<ly;++j)for(int i=0;i<lx;++i){const int n=((e*lz+k)*ly+j)*lx+i;
        for(int c=0;c<5;++c)plane[c]+=res[c][n];}
      for(int i=0;i<lx;++i){T cum[5]={0,0,0,0,0};
        for(int j=0;j<ly-1;++j){T line[5]={0,0,0,0,0};
          for(int ii=0;ii<lx;++ii){const int n=((e*lz+k)*ly+j)*lx+ii;
            for(int c=0;c<5;++c)line[c]+=res[c][n];}
          for(int c=0;c<5;++c){cum[c]+=wx[i]/sx*(line[c]-wy[j]/sy*plane[c]);
            correction[5*edge+c]=cum[c];} ++edge;}
      }
    }
    edge=ybase+lx*(ly-1)*lz;
    if(lz>1)for(int j=0;j<ly;++j)for(int i=0;i<lx;++i){T cum[5]={0,0,0,0,0};
      for(int k=0;k<lz-1;++k){T plane[5]={0,0,0,0,0};
        for(int jj=0;jj<ly;++jj)for(int ii=0;ii<lx;++ii){const int n=((e*lz+k)*ly+jj)*lx+ii;
          for(int c=0;c<5;++c)plane[c]+=res[c][n];}
        for(int c=0;c<5;++c){cum[c]+=wx[i]/sx*wy[j]/sy*
            (plane[c]-wz[k]/sz*element_sum[c]); correction[5*edge+c]=cum[c];}
        ++edge;}
    }
  }
}

template<typename T>
__device__ void line_node_residual(int direction,int pos,int i,int j,int k,int e,
  int lx,int ly,int lz,int edge_first,const T *rho,const T *mx,const T *my,
  const T *mz,const T *energy,const T *coef,const T *B,const T *jacinv,
  const T *mr0,const T *mr1,const T *mr2,const T *D,T gamma,T low[5],T high[5]) {
  const int length=direction==0?lx:(direction==1?ly:lz);
  int ci=i,cj=j,ck=k; if(direction==0)ci=pos; else if(direction==1)cj=pos; else ck=pos;
  const int node=((e*lz+ck)*ly+cj)*lx+ci;
  for(int c=0;c<5;++c){low[c]=T(0);high[c]=T(0);}
  if(pos<length-1){int ni=ci,nj=cj,nk=ck;if(direction==0)++ni;else if(direction==1)++nj;else ++nk;
    int nn=((e*lz+nk)*ly+nj)*lx+ni;T s[5],x[5],y[5],z[5];load_state(rho,mx,my,mz,energy,nn,s);fluxes(s,gamma,x,y,z);
    int ed=edge_first+pos;for(int c=0;c<5;++c)low[c]+=coef[3*ed]*x[c]+coef[3*ed+1]*y[c]+coef[3*ed+2]*z[c];}
  if(pos>0){int ni=ci,nj=cj,nk=ck;if(direction==0)--ni;else if(direction==1)--nj;else --nk;
    int nn=((e*lz+nk)*ly+nj)*lx+ni;T s[5],x[5],y[5],z[5];load_state(rho,mx,my,mz,energy,nn,s);fluxes(s,gamma,x,y,z);
    int ed=edge_first+pos-1;for(int c=0;c<5;++c)low[c]-=coef[3*ed]*x[c]+coef[3*ed+1]*y[c]+coef[3*ed+2]*z[c];}
  if(pos==0||pos==length-1){T s[5],x[5],y[5],z[5];load_state(rho,mx,my,mz,energy,node,s);fluxes(s,gamma,x,y,z);
    int ed=pos==0?edge_first:edge_first+length-2;const T sign=pos==0?T(-1):T(1);
    for(int c=0;c<5;++c)low[c]+=sign*(coef[3*ed]*x[c]+coef[3*ed+1]*y[c]+coef[3*ed+2]*z[c]);}
  for(int l=0;l<length;++l){int ni=i,nj=j,nk=k;if(direction==0)ni=l;else if(direction==1)nj=l;else nk=l;
    int nn=((e*lz+nk)*ly+nj)*lx+ni;T s[5],x[5],y[5],z[5];load_state(rho,mx,my,mz,energy,nn,s);fluxes(s,gamma,x,y,z);
    const T d=D[pos+l*length];for(int c=0;c<5;++c)high[c]+=d*(mr0[node]*x[c]+mr1[node]*y[c]+mr2[node]*z[c]);}
  for(int c=0;c<5;++c)high[c]*=B[node]*jacinv[node];
}

template<typename T>
__global__ void reconstruct_direction_kernel(int direction,const T *rho,
  const T *mx,const T *my,const T *mz,const T *energy,const T *coef,const T *B,
  const T *jacinv,const T *mr0,const T *mr1,const T *mr2,const T *D,
  T *correction,T *error,T gamma,int lx,int ly,int lz,int edges_per_element,
  int line_count) {
  for(int line=blockIdx.x*blockDim.x+threadIdx.x;line<line_count;
      line+=blockDim.x*gridDim.x){int q=line,i=0,j=0,k=0,e=0,first=0,length=0;
    const int xcount=(lx-1)*ly*lz,ycount=lx*(ly-1)*lz;
    if(direction==0){j=q%ly;q/=ly;k=q%lz;e=q/lz;length=lx;first=e*edges_per_element+(k*ly+j)*(lx-1);}
    else if(direction==1){i=q%lx;q/=lx;k=q%lz;e=q/lz;length=ly;first=e*edges_per_element+xcount+(k*lx+i)*(ly-1);}
    else{i=q%lx;q/=lx;j=q%ly;e=q/ly;length=lz;first=e*edges_per_element+xcount+ycount+(j*lx+i)*(lz-1);}
    T cumulative[5]={0,0,0,0,0},sum[5]={0,0,0,0,0},scale[5]={0,0,0,0,0};
    for(int pos=0;pos<length;++pos){T low[5],high[5];line_node_residual(direction,pos,i,j,k,e,lx,ly,lz,first,
        rho,mx,my,mz,energy,coef,B,jacinv,mr0,mr1,mr2,D,gamma,low,high);
      for(int c=0;c<5;++c){const T value=low[c]-high[c];sum[c]+=value;scale[c]+=fabs(low[c])+fabs(high[c]);
        if(pos<length-1){cumulative[c]+=value;correction[5*(first+pos)+c]=cumulative[c];}}
    }
    for(int c=0;c<5;++c)atomic_max_real(error+c,fabs(sum[c])/rmax(T(1),scale[c]));
  }
}

template<typename T>
__global__ void floor_timestep_kernel(const T *rho,const T *mx,const T *my,
  const T *mz,const T *energy,const T *r0,const T *r1,const T *r2,const T *r3,
  const T *r4,T *limit,T floor,T upper,int n) {
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x){
    T s[5],r[5];load_state(rho,mx,my,mz,energy,i,s);load_state(r0,r1,r2,r3,r4,i,r);
    limit[i]=floor_timestep(s,r,floor,upper);}
}

template<typename T>
__global__ void low_update_kernel(const T *rho,const T *mx,const T *my,
  const T *mz,const T *energy,T *r0,T *r1,T *r2,T *r3,T *r4,T dt,
  int scalar_mode,int n) {T *out[5]={r0,r1,r2,r3,r4};const T *in[5]={rho,mx,my,mz,energy};
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x)
    for(int c=0;c<5;++c)out[c][i]=(scalar_mode&&c>0)?in[c][i]:in[c][i]-dt*out[c][i];
}

template<typename T>
__global__ void bounds_init_kernel(const T *rho,const T *mx,const T *my,
  const T *mz,const T *energy,T *lower,T *upper,T *entropy,T *work,T gamma,
  int use_entropy,int n) {for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x){
    lower[i]=upper[i]=rho[i];work[i]=T(0);if(use_entropy){T s[5];load_state(rho,mx,my,mz,energy,i,s);entropy[i]=specific_entropy(s,gamma);}}
}

template<typename T>
__global__ void bounds_edge_kernel(const T *rho,const T *mx,const T *my,
  const T *mz,const int *left,const int *right,const T *coef,const T *visc,
  const T *second_weight_left,const T *second_weight_right,T *lower,T *upper,
  const T *entropy_stage,T *entropy,T *work,int relax,int use_entropy,
  int nedge) {
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nedge;e+=blockDim.x*gridDim.x){
    const int a=left[e],b=right[e];const T ra=rho[a],rb=rho[b];
    const T fa=coef[3*e]*mx[a]+coef[3*e+1]*my[a]+coef[3*e+2]*mz[a];
    const T fb=coef[3*e]*mx[b]+coef[3*e+1]*my[b]+coef[3*e+2]*mz[b];
    T bar=T(0.5)*(ra+rb);if(visc[e]>real_tiny<T>())bar-=(fb-fa)/(T(2)*visc[e]);
    atomic_min_real(lower+a,rmin(rb,bar));atomic_min_real(lower+b,rmin(ra,bar));
    atomic_max_real(upper+a,rmax(rb,bar));atomic_max_real(upper+b,rmax(ra,bar));
    if(use_entropy){const T ea=entropy_stage[a],eb=entropy_stage[b];atomic_min_real(entropy+a,eb);atomic_min_real(entropy+b,ea);}
    if(relax){const T diff=ra-rb;atomicAdd(work+a,second_weight_left[e]*diff);atomicAdd(work+b,-second_weight_right[e]*diff);}
  }
}

template<typename T>
__global__ void relax_edge_kernel(const int *left,const int *right,
  const int *direction,const T *degree0,const T *degree1,const T *degree2,
  const T *first,T *second,int nedge){const T *degree[3]={degree0,degree1,degree2};
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nedge;e+=blockDim.x*gridDim.x){const int a=left[e],b=right[e],d=direction[e]-1;
    const T avg=T(0.5)*(first[a]+first[b]);atomicAdd(second+a,T(2)/degree[d][a]*avg);atomicAdd(second+b,T(2)/degree[d][b]*avg);}}

template<typename T>
__global__ void relax_finalize_kernel(T *lower,T *upper,const T *second,
  T factor,T nodal_mass,T volume,int dimensions,int n){const T rh=pow(nodal_mass/volume,T(1.5)/T(dimensions));
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x){const T d=second[i]/(T(2)*T(2*dimensions+1));
    const T relax=factor*fabs(d);lower[i]=rmax((T(1)-rh)*lower[i],lower[i]-relax);upper[i]+=relax;}}

template<typename T>
__global__ void entropy_relax_edge_kernel(const T *rho,const T *mx,
  const T *my,const T *mz,const T *energy,const int *left,const int *right,
  T *midpoint_max,T gamma,int nedge) {
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nedge;
      e+=blockDim.x*gridDim.x) {
    const int a=left[e],b=right[e];
    T sa[5],sb[5],midpoint[5];
    load_state(rho,mx,my,mz,energy,a,sa);
    load_state(rho,mx,my,mz,energy,b,sb);
    for(int c=0;c<5;++c) midpoint[c]=T(0.5)*(sa[c]+sb[c]);
    const T entropy=specific_entropy(midpoint,gamma);
    atomic_max_real(midpoint_max+a,entropy);
    atomic_max_real(midpoint_max+b,entropy);
  }
}

template<typename T>
__global__ void entropy_relax_finalize_kernel(T *lower,
  const T *midpoint_max,T factor,T nodal_mass,T volume,int dimensions,int n) {
  const T rh=pow(nodal_mass/volume,T(1.5)/T(dimensions));
  const T cutoff=log(T(1)+rh/factor);
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;
      i+=blockDim.x*gridDim.x) {
    const T difference=rmax(T(0),midpoint_max[i]-lower[i]);
    const T variation=exp(rmin(difference,cutoff))-T(1);
    const T relaxation=rmin(rh,
                            factor*rmax(T(0),variation));
    lower[i]+=log(T(1)-relaxation);
  }
}

template<typename T>
__global__ void blend_kernel(const T *rho,const T *mx,const T *my,const T *mz,
  const T *energy,const int *left,const int *right,const T *visc,const T *ev,
  T *correction,T dt,int has_ev,int low_only,int scalar_mode,int nedge){const T *u[5]={rho,mx,my,mz,energy};
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nedge;e+=blockDim.x*gridDim.x){const int a=left[e],b=right[e];T high=low_only?T(0):T(1);
    if(!low_only&&has_ev)high=T(1)-rmin(T(1),rmax(T(0),rmax(ev[a],ev[b])));
    for(int c=0;c<5;++c)correction[5*e+c]=(scalar_mode&&c>0)?T(0):dt*high*(correction[5*e+c]+visc[e]*(u[c][a]-u[c][b]));}}

template<typename T>
__global__ void limiter_kernel(const T *q0,const T *q1,const T *q2,
  const T *q3,const T *q4,const T *lower,const T *upper,const T *entropy,
  const T *mass,const T *degree0,const T *degree1,const T *degree2,
  const int *left,const int *right,const int *direction,T *correction,
  T *edge_limit,T *limited,T *density_flag,T *energy_flag,T *entropy_flag,
  T gamma,T floor,int enforce_energy,int enforce_entropy,int check_base,
  int diagnostics_level,int dimensions,int nedge){const T *q[5]={q0,q1,q2,q3,q4};const T *degree[3]={degree0,degree1,degree2};
  const T eps=real_epsilon<T>();
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nedge;e+=blockDim.x*gridDim.x){const int a=left[e],b=right[e],d=direction[e]-1;
    T lb[5],rb[5],lc[5],rc[5];T scale=T(1),size=T(0);for(int c=0;c<5;++c){lb[c]=q[c][a];rb[c]=q[c][b];
      lc[c]=correction[5*e+c]*T(dimensions)*degree[d][a]/mass[a];rc[c]=-correction[5*e+c]*T(dimensions)*degree[d][b]/mass[b];
      scale=rmax(scale,rmax(fabs(lb[c]),fabs(rb[c])));size=rmax(size,rmax(fabs(lc[c]),fabs(rc[c])));}
    bool dl=false,el=false,sl=false;T limit=T(1);
    if(size<=T(512)*eps*scale){for(int c=0;c<5;++c)correction[5*e+c]=T(0);}
    else {bool dl1,el1,sl1,dl2,el2,sl2;const T sent=enforce_entropy?entropy[a]:-real_huge<T>();
      const T tent=enforce_entropy?entropy[b]:-real_huge<T>();
      T l1=limit_endpoint(lb,lc,lower[a],upper[a],sent,gamma,floor,enforce_energy,enforce_entropy,check_base,dl1,el1,sl1);
      T l2=limit_endpoint(rb,rc,lower[b],upper[b],tent,gamma,floor,enforce_energy,enforce_entropy,check_base,dl2,el2,sl2);
      limit=rmin(l1,l2);dl=dl1||dl2;el=el1||el2;sl=sl1||sl2;
      if(enforce_entropy&&limit>T(0)&&limit<T(1)){T lt[5],rt[5];for(int c=0;c<5;++c){lt[c]=lb[c]+limit*lc[c];rt[c]=rb[c]+limit*rc[c];}
        if(!entropy_is_admissible(lt,gamma,sent,floor)||!entropy_is_admissible(rt,gamma,tent,floor)){
          T lo=T(0),hi=limit;for(int it=0;it<64;++it){const T mid=lo+T(0.5)*(hi-lo);for(int c=0;c<5;++c){lt[c]=lb[c]+mid*lc[c];rt[c]=rb[c]+mid*rc[c];}
            if(entropy_is_admissible(lt,gamma,sent,floor)&&entropy_is_admissible(rt,gamma,tent,floor))lo=mid;else hi=mid;}limit=lo;sl=true;}}
      for(int c=0;c<5;++c)correction[5*e+c]*=limit;}
    if(diagnostics_level>0)edge_limit[e]=limit;
    if(diagnostics_level>1){limited[e]=limit<T(1)-T(32)*eps?T(1):T(0);
      density_flag[e]=dl?T(1):T(0);energy_flag[e]=el?T(1):T(0);
      entropy_flag[e]=sl?T(1):T(0);}}
}

template<typename T>
__global__ void incidence_kernel(const int *left,const int *right,
  const T *correction,T *r0,T *r1,T *r2,T *r3,T *r4,int nedge){T *r[5]={r0,r1,r2,r3,r4};
  for(int e=blockIdx.x*blockDim.x+threadIdx.x;e<nedge;e+=blockDim.x*gridDim.x){const int a=left[e],b=right[e];
    for(int c=0;c<5;++c){atomicAdd(r[c]+a,correction[5*e+c]);atomicAdd(r[c]+b,-correction[5*e+c]);}}}

template<typename T>
__global__ void correction_update_kernel(T *q0,T *q1,T *q2,T *q3,T *q4,
  const T *r0,const T *r1,const T *r2,const T *r3,const T *r4,const T *mass,int n){T *q[5]={q0,q1,q2,q3,q4};const T *r[5]={r0,r1,r2,r3,r4};
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x)for(int c=0;c<5;++c)q[c][i]+=r[c][i]/mass[i];}

template<typename T>
__global__ void validation_kernel(const T *q0,const T *q1,const T *q2,
  const T *q3,const T *q4,const T *lower,const T *upper,const T *entropy_lower,
  T *lower_v,T *upper_v,T *entropy_v,T *entropy_value,T gamma,int use_entropy,int n){
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x){T s[5];load_state(q0,q1,q2,q3,q4,i,s);
    lower_v[i]=rmax(T(0),lower[i]-q0[i]);upper_v[i]=rmax(T(0),q0[i]-upper[i]);
    const T ent=specific_entropy(s,gamma);entropy_value[i]=ent;
    entropy_v[i]=use_entropy?rmax(T(0),entropy_lower[i]-ent-entropy_tolerance(s,entropy_lower[i])):T(0);}}

enum diagnostic_summary_kind {
  full_summary = 0,
  validation_summary = 1,
  observation_summary = 2,
  limiter_status_summary = 3
};

template<typename T>
__global__ void diagnostic_summary_init_kernel(T *summary, int size,
                                                int kind) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  for (int i = 0; i < size; ++i) summary[i] = T(0);
  if (kind == full_summary) {
    summary[3] = real_huge<T>();
  } else if (kind == validation_summary) {
    for (int i = 3; i <= 6; ++i) summary[i] = real_huge<T>();
  } else if (kind == observation_summary) {
    for (int i = 0; i <= 2; ++i) summary[i] = real_huge<T>();
  } else if (kind == limiter_status_summary) {
    summary[0] = real_huge<T>();
    summary[1] = -real_huge<T>();
  }
}

enum full_diagnostic_partial_size {
  full_node_partial_size = 8,
  full_edge_partial_size = 10
};

template<typename T>
__global__ void full_diagnostics_node_partial_kernel(const T *rho,
  const T *mx, const T *my, const T *mz, const T *energy, T gamma,
  const T *entropy_fraction, int has_entropy, const T *r0, const T *r1,
  const T *r2, const T *r3, const T *r4, T *partial, int n) {
  extern __shared__ unsigned char storage[];
  T *shared = reinterpret_cast<T *>(storage);
  const T *residual[5] = {r0, r1, r2, r3, r4};
  const int tid = threadIdx.x;
  const int i = blockIdx.x * blockDim.x + tid;
  T value[full_node_partial_size] = {T(0), T(0), T(0), T(0),
                                     T(0), T(0), T(0), T(0)};
  if (i < n) {
    T state[5];
    load_state(rho, mx, my, mz, energy, i, state);
    const T pressure = (gamma - T(1)) * internal_energy(state);
    const T velocity = sqrt(mx[i]*mx[i] + my[i]*my[i] + mz[i]*mz[i]) /
      rho[i];
    value[0] = velocity + sqrt(gamma * pressure / rho[i]);
    if (has_entropy) {
      value[1] = entropy_fraction[i];
      value[2] = entropy_fraction[i];
    }
    for (int c = 0; c < 5; ++c) value[3+c] = residual[c][i];
  }
  for (int metric = 0; metric < full_node_partial_size; ++metric)
    shared[metric*blockDim.x+tid] = value[metric];
  __syncthreads();
  for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
    if (tid < stride) {
      shared[tid] = rmax(shared[tid], shared[tid+stride]);
      shared[blockDim.x+tid] = rmax(shared[blockDim.x+tid],
                                    shared[blockDim.x+tid+stride]);
      for (int metric = 2; metric < full_node_partial_size; ++metric)
        shared[metric*blockDim.x+tid] +=
          shared[metric*blockDim.x+tid+stride];
    }
    __syncthreads();
  }
  if (tid < full_node_partial_size)
    partial[blockIdx.x*full_node_partial_size+tid] =
      shared[tid*blockDim.x];
}

template<typename T>
__global__ void full_diagnostics_edge_partial_kernel(const T *edge_limit,
  const T *limited, const T *density_flag, const T *energy_flag,
  const T *entropy_flag, const T *correction, T *partial, int nedge) {
  extern __shared__ unsigned char storage[];
  T *shared = reinterpret_cast<T *>(storage);
  const int tid = threadIdx.x;
  const int i = blockIdx.x * blockDim.x + tid;
  const T huge = real_huge<T>();
  const T eps = real_epsilon<T>();
  T value[full_edge_partial_size] = {huge, -huge, T(0), T(0), T(0),
                                     T(0), T(0), T(0), T(0), T(0)};
  if (i < nedge) {
    const T limit = edge_limit[i];
    value[0] = limit;
    value[1] = limit;
    value[2] = limit;
    value[3] = limited[i];
    value[4] = density_flag[i];
    value[5] = energy_flag[i];
    value[6] = entropy_flag[i];
    value[7] = (!isfinite(limit) || limit < -T(32)*eps ||
                limit > T(1) + T(32)*eps) ? T(1) : T(0);
    for (int c = 0; c < 5; ++c) {
      const T flux = correction[5*i+c];
      value[8] = rmax(value[8], fabs(flux));
      value[9] += flux*flux;
    }
  }
  for (int metric = 0; metric < full_edge_partial_size; ++metric)
    shared[metric*blockDim.x+tid] = value[metric];
  __syncthreads();
  for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
    if (tid < stride) {
      shared[tid] = rmin(shared[tid], shared[tid+stride]);
      shared[blockDim.x+tid] = rmax(shared[blockDim.x+tid],
                                    shared[blockDim.x+tid+stride]);
      for (int metric = 2; metric < 8; ++metric)
        shared[metric*blockDim.x+tid] +=
          shared[metric*blockDim.x+tid+stride];
      shared[8*blockDim.x+tid] = rmax(shared[8*blockDim.x+tid],
                                      shared[8*blockDim.x+tid+stride]);
      shared[9*blockDim.x+tid] += shared[9*blockDim.x+tid+stride];
    }
    __syncthreads();
  }
  if (tid < full_edge_partial_size)
    partial[blockIdx.x*full_edge_partial_size+tid] =
      shared[tid*blockDim.x];
}

template<typename T>
__global__ void full_diagnostics_finalize_kernel(const T *node_partial,
  const T *edge_partial, const T *directional_error, T *summary,
  int node_blocks, int edge_blocks) {
  const int metric = threadIdx.x;
  if (blockIdx.x != 0 || metric >= 23) return;
  T result = T(0);
  if (metric == 0 || metric == 1) {
    for (int block = 0; block < node_blocks; ++block)
      result = rmax(result,
        node_partial[block*full_node_partial_size+metric]);
  } else if (metric == 2) {
    for (int block = 0; block < node_blocks; ++block)
      result += node_partial[block*full_node_partial_size+2];
  } else if (metric >= 10 && metric < 15) {
    const int source = metric - 7;
    for (int block = 0; block < node_blocks; ++block)
      result += node_partial[block*full_node_partial_size+source];
  } else if (metric >= 15 && metric < 20) {
    result = directional_error[metric-15];
  } else {
    int source = -1;
    if (metric >= 3 && metric <= 9) source = metric - 3;
    else if (metric == 20) source = 8;
    else if (metric == 21) source = 9;
    else if (metric == 22) source = 7;
    if (metric == 3) result = real_huge<T>();
    else if (metric == 4) result = -real_huge<T>();
    for (int block = 0; block < edge_blocks; ++block) {
      const T value = edge_partial[block*full_edge_partial_size+source];
      if (metric == 3) result = rmin(result, value);
      else if (metric == 4 || metric == 20) result = rmax(result, value);
      else result += value;
    }
  }
  summary[metric] = result;
}

template<typename T>
__global__ void validation_summary_kernel(const T *q0, const T *q1,
  const T *q2, const T *q3, const T *q4, const T *lower, const T *upper,
  const T *entropy_lower, T *summary, T gamma, int use_entropy, int n) {
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += blockDim.x * gridDim.x) {
    T state[5];
    load_state(q0, q1, q2, q3, q4, i, state);
    const T internal = internal_energy(state);
    const T pressure = (gamma - T(1)) * internal;
    const T entropy = specific_entropy(state, gamma);
    const T lower_violation = rmax(T(0), lower[i] - q0[i]);
    const T upper_violation = rmax(T(0), q0[i] - upper[i]);
    const T entropy_violation = use_entropy ?
      rmax(T(0), entropy_lower[i] - entropy -
           entropy_tolerance(state, entropy_lower[i])) : T(0);
    atomic_max_real(summary, lower_violation);
    atomic_max_real(summary + 1, upper_violation);
    atomic_max_real(summary + 2, entropy_violation);
    atomic_min_real(summary + 3, q0[i]);
    atomic_min_real(summary + 4, internal);
    atomic_min_real(summary + 5, pressure);
    atomic_min_real(summary + 6, entropy);
    atomic_max_real(summary + 7, rmax(q0[i], upper[i]));
    bool invalid = q0[i] <= T(0) || !isfinite(internal) || internal <= T(0) ||
      !isfinite(pressure) || pressure <= T(0) || !isfinite(entropy);
    for (int c = 0; c < 5; ++c) invalid = invalid || !isfinite(state[c]);
    if (invalid) atomicAdd(summary + 8, T(1));
  }
}

template<typename T>
__global__ void observation_summary_kernel(const T *rho, const T *internal,
  const T *pressure, const T *u, const T *v, const T *w, const T *sound,
  T *summary, int n) {
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += blockDim.x * gridDim.x) {
    atomic_min_real(summary, rho[i]);
    atomic_min_real(summary + 1, internal[i]);
    atomic_min_real(summary + 2, pressure[i]);
    const T wave = sqrt(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]) + sound[i];
    atomic_max_real(summary + 3, wave);
  }
}

template<typename T>
__global__ void limiter_status_kernel(const T *edge_limit, T *summary,
                                      int nedge) {
  const T eps = real_epsilon<T>();
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nedge;
       i += blockDim.x * gridDim.x) {
    const T limit = edge_limit[i];
    atomic_min_real(summary, limit);
    atomic_max_real(summary + 1, limit);
    if (!isfinite(limit) || limit < -T(32)*eps ||
        limit > T(1) + T(32)*eps)
      atomicAdd(summary + 2, T(1));
  }
}

template<typename T>
__global__ void update_uvw_kernel_idp(T *u,T *v,T *w,const T *mx,const T *my,
  const T *mz,const T *rho,int n){for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x){u[i]=mx[i]/rho[i];v[i]=my[i]/rho[i];w[i]=mz[i]/rho[i];}}

template<typename T>
__global__ void update_momentum_pressure_kernel(T *mx,T *my,T *mz,T *p,T *kin,
  const T *u,const T *v,const T *w,const T *energy,const T *rho,T gamma,int n){
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=blockDim.x*gridDim.x){mx[i]=rho[i]*u[i];my[i]=rho[i]*v[i];mz[i]=rho[i]*w[i];
    kin[i]=T(0.5)*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]);p[i]=(gamma-T(1))*(energy[i]-kin[i]);}}

template<typename T>
__global__ void update_energy_kernel(T *energy, T *p, const T *kin, T gamma,
                                     T floor, int n) {
  const T inv = T(1) / (gamma - T(1));
  floor = rmax(floor, T(0));
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < n;
       i += blockDim.x*gridDim.x) {
    if (p[i]*inv <= floor) {
      p[i] = rmax(p[i], T(1e-12));
      if (p[i]*inv <= floor) {
        const T guard = T(32)*real_epsilon<T>()*
          rmax(T(1), rmax(fabs(kin[i]), floor));
        p[i] = (gamma - T(1))*(floor + guard);
      }
    }
    energy[i] = p[i]*inv + kin[i];
  }
}

inline dim3 blocks(int n) { return dim3((n+255)/256,1,1); }
inline dim3 threads() { return dim3(256,1,1); }

} // namespace

extern "C" {

void gpu_euler_idp_primitives(void *rho,void *mx,void *my,void *mz,void *energy,
  void *u,void *v,void *w,void *p,void *sound,void *internal,void *status,
  real *gamma,real *floor,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  primitives_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(real*)u,(real*)v,(real*)w,(real*)p,(real*)sound,(real*)internal,(real*)status,*gamma,*floor,*n);idp_check_launch();}

void gpu_euler_idp_graph_viscosity(void *rho,void *mx,void *my,void *mz,
  void *energy,void *user_speed,int *has_user,void *left,void *right,void *coef,
  void *edge_visc,void *vis_sum,void *edge_speed,real *gamma,int *n,int *nedge){idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  idp_zero(vis_sum,sizeof(real)*(*n),s);graph_viscosity_kernel<real><<<blocks(*nedge),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(real*)user_speed,*has_user,(int*)left,(int*)right,(real*)coef,(real*)edge_visc,(real*)vis_sum,(real*)edge_speed,*gamma,*nedge);idp_check_launch();}

void gpu_euler_idp_graph_rate(void *vis,void *mass,void *rate,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;graph_rate_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)vis,(real*)mass,(real*)rate,*n);idp_check_launch();}
void gpu_euler_idp_nodal_wave(void *u,void *v,void *w,void *sound,void *wave,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;nodal_wave_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)u,(real*)v,(real*)w,(real*)sound,(real*)wave,*n);idp_check_launch();}

void gpu_euler_idp_low_residual(void *rho,void *mx,void *my,void *mz,void *energy,
  void *diag,void *left,void *right,void *coef,void *visc,void *r0,void *r1,
  void *r2,void *r3,void *r4,real *gamma,int *periodic,int *n,int *nedge){idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  low_init_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(real*)diag,(real*)r0,(real*)r1,(real*)r2,(real*)r3,(real*)r4,*gamma,*periodic,*n);
  low_edge_kernel<real><<<blocks(*nedge),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(int*)left,(int*)right,(real*)coef,(real*)visc,(real*)r0,(real*)r1,(real*)r2,(real*)r3,(real*)r4,*gamma,*periodic,*nedge);idp_check_launch();}

void gpu_euler_idp_scale_residual(void *r0,void *r1,void *r2,void *r3,void *r4,
  void *binv,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;scale_residual_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)r0,(real*)r1,(real*)r2,(real*)r3,(real*)r4,(real*)binv,*n);idp_check_launch();}

void gpu_euler_idp_reconstruct(
    void *rho, void *mx, void *my, void *mz, void *energy, void *diag,
    void *left, void *right, void *coef, void *const residual[5],
    void *correction, void *error, void *B, void *jacinv, void *const dr[3],
    void *const ds[3], void *const dt[3], void *const derivative[3],
    void *const weight[3], real *gamma, int *lx, int *ly, int *lz, int *nelv,
    int *edges_per_element, int *affine, int *n, int *nedge) {
  idp_stream_t s = (idp_stream_t)glb_cmd_queue;
  idp_zero(correction, sizeof(real) * 5 * (*nedge), s);
  idp_zero(error, sizeof(real) * 5, s);
  if (*affine) {
    const int nx = (*nelv) * (*ly) * (*lz), ny = (*nelv) * (*lx) * (*lz),
              nz = (*nelv) * (*lx) * (*ly);
    reconstruct_direction_kernel<real><<<blocks(nx), threads(), 0, s>>>(
        0, (real *)rho, (real *)mx, (real *)my, (real *)mz, (real *)energy,
        (real *)coef, (real *)B, (real *)jacinv, (real *)dr[0], (real *)dr[1],
        (real *)dr[2], (real *)derivative[0], (real *)correction, (real *)error,
        *gamma, *lx, *ly, *lz, *edges_per_element, nx);
    reconstruct_direction_kernel<real><<<blocks(ny), threads(), 0, s>>>(
        1, (real *)rho, (real *)mx, (real *)my, (real *)mz, (real *)energy,
        (real *)coef, (real *)B, (real *)jacinv, (real *)ds[0], (real *)ds[1],
        (real *)ds[2], (real *)derivative[1], (real *)correction, (real *)error,
        *gamma, *lx, *ly, *lz, *edges_per_element, ny);
    if (*lz > 1)
      reconstruct_direction_kernel<real><<<blocks(nz), threads(), 0, s>>>(
          2, (real *)rho, (real *)mx, (real *)my, (real *)mz, (real *)energy,
          (real *)coef, (real *)B, (real *)jacinv, (real *)dt[0], (real *)dt[1],
          (real *)dt[2], (real *)derivative[2], (real *)correction,
          (real *)error, *gamma, *lx, *ly, *lz, *edges_per_element, nz);
  } else {
    correction_residual_init_kernel<real><<<blocks(*n), threads(), 0, s>>>(
        (real *)rho, (real *)mx, (real *)my, (real *)mz, (real *)energy,
        (real *)diag, (real *)dr[0], (real *)dr[1], (real *)dr[2],
        (real *)ds[0], (real *)ds[1], (real *)ds[2], (real *)dt[0],
        (real *)dt[1], (real *)dt[2], (real *)derivative[0],
        (real *)derivative[1], (real *)derivative[2], (real *)weight[0],
        (real *)weight[1], (real *)weight[2], (real *)residual[0],
        (real *)residual[1], (real *)residual[2], (real *)residual[3],
        (real *)residual[4], *gamma, *lx, *ly, *lz, *n);
    correction_residual_edge_kernel<real><<<blocks(*nedge), threads(), 0, s>>>(
        (real *)rho, (real *)mx, (real *)my, (real *)mz, (real *)energy,
        (int *)left, (int *)right, (real *)coef, (real *)residual[0],
        (real *)residual[1], (real *)residual[2], (real *)residual[3],
        (real *)residual[4], *gamma, *nedge);
    reconstruct_general_kernel<real><<<blocks(*nelv), threads(), 0, s>>>(
        (real *)residual[0], (real *)residual[1], (real *)residual[2],
        (real *)residual[3], (real *)residual[4], (real *)weight[0],
        (real *)weight[1], (real *)weight[2], (real *)correction, (real *)error,
        *lx, *ly, *lz, *edges_per_element, *nelv);
  }
  idp_check_launch();
}

void gpu_euler_idp_floor_timestep(void *rho,void *mx,void *my,void *mz,
  void *energy,void *r0,void *r1,void *r2,void *r3,void *r4,void *limit,
  real *floor,real *upper,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;floor_timestep_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(real*)r0,(real*)r1,(real*)r2,(real*)r3,(real*)r4,(real*)limit,*floor,*upper,*n);idp_check_launch();}
void gpu_euler_idp_low_update(void *rho,void *mx,void *my,void *mz,void *energy,
  void *r0,void *r1,void *r2,void *r3,void *r4,real *dt,int *scalar_mode,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;low_update_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(real*)r0,(real*)r1,(real*)r2,(real*)r3,(real*)r4,*dt,*scalar_mode,*n);idp_check_launch();}

void gpu_euler_idp_bounds_init(void *rho,void *mx,void *my,void *mz,void *energy,
  void *lower,void *upper,void *entropy,void *work,real *gamma,int *use_entropy,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;bounds_init_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(real*)lower,(real*)upper,(real*)entropy,(real*)work,*gamma,*use_entropy,*n);idp_check_launch();}
void gpu_euler_idp_bounds_edges(void *rho,void *mx,void *my,void *mz,void *left,
  void *right,void *coef,void *visc,void *second_weight_left,
  void *second_weight_right,void *lower,void *upper,void *entropy_stage,
  void *entropy,void *work,int *relax,int *use_entropy,int *nedge){idp_stream_t s=(idp_stream_t)glb_cmd_queue;bounds_edge_kernel<real><<<blocks(*nedge),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(int*)left,(int*)right,(real*)coef,(real*)visc,(real*)second_weight_left,(real*)second_weight_right,(real*)lower,(real*)upper,(real*)entropy_stage,(real*)entropy,(real*)work,*relax,*use_entropy,*nedge);idp_check_launch();}
void gpu_euler_idp_relax_edges(void *left,void *right,void *direction,
  void *degree0,void *degree1,void *degree2,void *first,void *second,int *n,
  int *nedge){idp_stream_t s=(idp_stream_t)glb_cmd_queue;idp_zero(second,sizeof(real)*(*n),s);relax_edge_kernel<real><<<blocks(*nedge),threads(),0,s>>>((int*)left,(int*)right,(int*)direction,(real*)degree0,(real*)degree1,(real*)degree2,(real*)first,(real*)second,*nedge);idp_check_launch();}
void gpu_euler_idp_relax_finalize(void *lower,void *upper,void *second,
  real *factor,real *nodal_mass,real *volume,int *dimensions,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;relax_finalize_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)lower,(real*)upper,(real*)second,*factor,*nodal_mass,*volume,*dimensions,*n);idp_check_launch();}

void gpu_euler_idp_entropy_relax_edges(void *rho,void *mx,void *my,void *mz,
  void *energy,void *left,void *right,void *midpoint_max,real *gamma,
  int *nedge){idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  entropy_relax_edge_kernel<real><<<blocks(*nedge),threads(),0,s>>>(
    (real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(int*)left,
    (int*)right,(real*)midpoint_max,*gamma,*nedge);
  idp_check_launch();}

void gpu_euler_idp_entropy_relax_finalize(void *lower,void *midpoint_max,
  real *factor,real *nodal_mass,real *volume,int *dimensions,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  entropy_relax_finalize_kernel<real><<<blocks(*n),threads(),0,s>>>(
    (real*)lower,(real*)midpoint_max,*factor,*nodal_mass,*volume,*dimensions,*n);
  idp_check_launch();}

void gpu_euler_idp_blend(void *rho,void *mx,void *my,void *mz,void *energy,
  void *left,void *right,void *visc,void *ev,void *correction,real *dt,
  int *has_ev,int *low_only,int *scalar_mode,int *nedge){idp_stream_t s=(idp_stream_t)glb_cmd_queue;blend_kernel<real><<<blocks(*nedge),threads(),0,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,(int*)left,(int*)right,(real*)visc,(real*)ev,(real*)correction,*dt,*has_ev,*low_only,*scalar_mode,*nedge);idp_check_launch();}
void gpu_euler_idp_limiter(void *q0,void *q1,void *q2,void *q3,void *q4,
  void *lower,void *upper,void *entropy,void *mass,void *degree0,void *degree1,
  void *degree2,void *left,void *right,void *direction,void *correction,
  void *edge_limit,void *limited,void *density_flag,void *energy_flag,
  void *entropy_flag,real *gamma,real *floor,int *enforce_energy,
  int *enforce_entropy,int *check_base,int *diagnostics_level,int *dimensions,
  int *nedge){idp_stream_t s=(idp_stream_t)glb_cmd_queue;limiter_kernel<real><<<blocks(*nedge),threads(),0,s>>>((real*)q0,(real*)q1,(real*)q2,(real*)q3,(real*)q4,(real*)lower,(real*)upper,(real*)entropy,(real*)mass,(real*)degree0,(real*)degree1,(real*)degree2,(int*)left,(int*)right,(int*)direction,(real*)correction,(real*)edge_limit,(real*)limited,(real*)density_flag,(real*)energy_flag,(real*)entropy_flag,*gamma,*floor,*enforce_energy,*enforce_entropy,*check_base,*diagnostics_level,*dimensions,*nedge);idp_check_launch();}
void gpu_euler_idp_incidence(void *left,void *right,void *correction,void *q0,
  void *q1,void *q2,void *q3,void *q4,void *r0,void *r1,void *r2,void *r3,
  void *r4,int *n,int *nedge){idp_stream_t s=(idp_stream_t)glb_cmd_queue;void *r[5]={r0,r1,r2,r3,r4};for(int c=0;c<5;++c)idp_zero(r[c],sizeof(real)*(*n),s);incidence_kernel<real><<<blocks(*nedge),threads(),0,s>>>((int*)left,(int*)right,(real*)correction,(real*)r0,(real*)r1,(real*)r2,(real*)r3,(real*)r4,*nedge);idp_check_launch();}
void gpu_euler_idp_correction_update(void *q0,void *q1,void *q2,void *q3,
  void *q4,void *r0,void *r1,void *r2,void *r3,void *r4,void *mass,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;correction_update_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)q0,(real*)q1,(real*)q2,(real*)q3,(real*)q4,(real*)r0,(real*)r1,(real*)r2,(real*)r3,(real*)r4,(real*)mass,*n);idp_check_launch();}

void gpu_euler_idp_validate(void *q0,void *q1,void *q2,void *q3,void *q4,
  void *lower,void *upper,void *entropy_lower,void *lower_v,void *upper_v,
  void *entropy_v,void *entropy_value,real *gamma,int *use_entropy,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;validation_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)q0,(real*)q1,(real*)q2,(real*)q3,(real*)q4,(real*)lower,(real*)upper,(real*)entropy_lower,(real*)lower_v,(real*)upper_v,(real*)entropy_v,(real*)entropy_value,*gamma,*use_entropy,*n);idp_check_launch();}

void gpu_euler_idp_full_diagnostics(void *rho,void *mx,void *my,void *mz,
  void *energy,real *gamma,void *entropy_fraction,int *has_entropy,
  void *edge_limit,void *limited,
  void *density_flag,void *energy_flag,void *entropy_flag,void *correction,
  void *r0,void *r1,void *r2,void *r3,void *r4,void *directional_error,
  void *summary,void *partial,int *n,int *nedge){
  idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  const int node_blocks=blocks(*n).x,edge_blocks=blocks(*nedge).x;
  real *node_partial=(real*)partial;
  real *edge_partial=node_partial+full_node_partial_size*node_blocks;
  const size_t node_shared=sizeof(real)*full_node_partial_size*threads().x;
  const size_t edge_shared=sizeof(real)*full_edge_partial_size*threads().x;
  full_diagnostics_node_partial_kernel<real><<<node_blocks,threads(),
    node_shared,s>>>((real*)rho,(real*)mx,(real*)my,(real*)mz,(real*)energy,
    *gamma,(real*)entropy_fraction,*has_entropy,(real*)r0,(real*)r1,
    (real*)r2,(real*)r3,(real*)r4,node_partial,*n);
  full_diagnostics_edge_partial_kernel<real><<<edge_blocks,threads(),
    edge_shared,s>>>((real*)edge_limit,(real*)limited,(real*)density_flag,
    (real*)energy_flag,(real*)entropy_flag,(real*)correction,edge_partial,
    *nedge);
  full_diagnostics_finalize_kernel<real><<<1,32,0,s>>>(node_partial,
    edge_partial,(real*)directional_error,(real*)summary,node_blocks,
    edge_blocks);
  idp_check_launch();}

void gpu_euler_idp_validation_summary(void *q0,void *q1,void *q2,void *q3,
  void *q4,void *lower,void *upper,void *entropy_lower,void *summary,
  real *gamma,int *use_entropy,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  diagnostic_summary_init_kernel<real><<<1,1,0,s>>>((real*)summary,9,validation_summary);
  validation_summary_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)q0,
    (real*)q1,(real*)q2,(real*)q3,(real*)q4,(real*)lower,(real*)upper,
    (real*)entropy_lower,(real*)summary,*gamma,*use_entropy,*n);
  idp_check_launch();}

void gpu_euler_idp_observation_summary(void *rho,void *internal,
  void *pressure,void *u,void *v,void *w,void *sound,void *summary,int *n){
  idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  diagnostic_summary_init_kernel<real><<<1,1,0,s>>>((real*)summary,4,observation_summary);
  observation_summary_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)rho,
    (real*)internal,(real*)pressure,(real*)u,(real*)v,(real*)w,(real*)sound,
    (real*)summary,*n);idp_check_launch();}

void gpu_euler_idp_limiter_status(void *edge_limit,void *summary,int *nedge){
  idp_stream_t s=(idp_stream_t)glb_cmd_queue;
  diagnostic_summary_init_kernel<real><<<1,1,0,s>>>((real*)summary,3,limiter_status_summary);
  limiter_status_kernel<real><<<blocks(*nedge),threads(),0,s>>>(
    (real*)edge_limit,(real*)summary,*nedge);idp_check_launch();}

void gpu_euler_idp_update_uvw(void *u,void *v,void *w,void *mx,void *my,
  void *mz,void *rho,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;update_uvw_kernel_idp<real><<<blocks(*n),threads(),0,s>>>((real*)u,(real*)v,(real*)w,(real*)mx,(real*)my,(real*)mz,(real*)rho,*n);idp_check_launch();}
void gpu_euler_idp_update_momentum_pressure(void *mx,void *my,void *mz,void *p,
  void *kin,void *u,void *v,void *w,void *energy,void *rho,real *gamma,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;update_momentum_pressure_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)mx,(real*)my,(real*)mz,(real*)p,(real*)kin,(real*)u,(real*)v,(real*)w,(real*)energy,(real*)rho,*gamma,*n);idp_check_launch();}
void gpu_euler_idp_update_energy(void *energy,void *p,void *kin,real *gamma,
  real *floor,int *n){idp_stream_t s=(idp_stream_t)glb_cmd_queue;update_energy_kernel<real><<<blocks(*n),threads(),0,s>>>((real*)energy,(real*)p,(real*)kin,*gamma,*floor,*n);idp_check_launch();}

}
