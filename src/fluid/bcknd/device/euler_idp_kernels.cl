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

/* The solver uses one numerical implementation across all GPU backends. */
#if defined(__METAL_VERSION__)
#include <metal_stdlib>
using namespace metal;
typedef float real;
#define IDP_GLOBAL device
#define IDP_PRIVATE thread
#define IDP_SHARED threadgroup
#define IDP_VALUE(type, name) constant type &name
#define IDP_GRID                                                                               \
  , uint IDP_GROUP [[threadgroup_position_in_grid]],                                           \
      uint IDP_WIDTH [[threads_per_threadgroup]], uint IDP_GROUPS [[threadgroups_per_grid]],   \
      uint IDP_LANE [[thread_index_in_threadgroup]]
#define IDP_BARRIER() threadgroup_barrier(mem_flags::mem_threadgroup)
#define IDP_KERNEL kernel
#define IDP_INLINE static inline
#define IDP_EPS FLT_EPSILON
#define IDP_TINY FLT_MIN
#define IDP_HUGE FLT_MAX
#elif defined(__OPENCL_C_VERSION__)
#pragma OPENCL FP_CONTRACT OFF
/* cltostring substitutes the configured floating-point type. */
#define IDP_IS_DOUBLE_float 0
#define IDP_IS_DOUBLE_double 1
#define IDP_IS_DOUBLE(type) IDP_IS_DOUBLE_##type
#if IDP_IS_DOUBLE(real)
#ifndef cl_khr_int64_base_atomics
#error "Double-precision Euler IDP requires 64-bit compare-and-swap"
#endif
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
typedef ulong idp_bits;
#define IDP_BITS as_ulong
#define IDP_FLOAT as_double
#define IDP_EPS DBL_EPSILON
#define IDP_TINY DBL_MIN
#define IDP_HUGE DBL_MAX
#else
typedef uint idp_bits;
#define IDP_BITS as_uint
#define IDP_FLOAT as_float
#define IDP_EPS FLT_EPSILON
#define IDP_TINY FLT_MIN
#define IDP_HUGE FLT_MAX
#endif
#define IDP_GLOBAL __global
#define IDP_PRIVATE __private
#define IDP_SHARED __local
#define IDP_VALUE(type, name) type name
#define IDP_GRID
#define IDP_GROUP ((int)get_group_id(0))
#define IDP_WIDTH ((int)get_local_size(0))
#define IDP_GROUPS ((int)get_num_groups(0))
#define IDP_LANE ((int)get_local_id(0))
#define IDP_BARRIER() barrier(CLK_LOCAL_MEM_FENCE)
#define IDP_KERNEL __kernel
#define IDP_INLINE static inline
#endif

IDP_INLINE real rmin(real a, real b) { return a < b ? a : b; }
IDP_INLINE real rmax(real a, real b) { return a > b ? a : b; }

#if defined(__METAL_VERSION__) || defined(__OPENCL_C_VERSION__)
/* Several edges can contribute to the same GLL node. Compare bit patterns
   so the retry also terminates for NaNs and signed zero. */
enum idp_atomic_op { idp_sum, idp_min, idp_max };
IDP_INLINE void idp_atomic(IDP_GLOBAL real *p, real value, enum idp_atomic_op op) {
#ifdef __METAL_VERSION__
  device atomic_uint *address = (device atomic_uint *)p;
  uint old = atomic_load_explicit(address, memory_order_relaxed);
#else
  volatile __global idp_bits *address = (volatile __global idp_bits *)p;
  idp_bits old = atom_cmpxchg(address, (idp_bits)0, (idp_bits)0);
#endif
  for (;;) {
#ifdef __METAL_VERSION__
    real current = as_type<real>(old);
#else
    real current = IDP_FLOAT(old);
#endif
    if ((op == idp_min && !(value < current)) || (op == idp_max && !(value > current)))
      return;
    real next = op == idp_sum ? current + value : value;
#ifdef __METAL_VERSION__
    if (atomic_compare_exchange_weak_explicit(address, &old, as_type<uint>(next),
                                              memory_order_relaxed, memory_order_relaxed))
      return;
#else
    idp_bits expected = old;
    old = atom_cmpxchg(address, expected, IDP_BITS(next));
    if (old == expected)
      return;
#endif
  }
}
#define idp_add(p, v) idp_atomic(p, v, idp_sum)
#define idp_minimum(p, v) idp_atomic(p, v, idp_min)
#define idp_maximum(p, v) idp_atomic(p, v, idp_max)

IDP_KERNEL void idp_zero_kernel(IDP_GLOBAL real *out, IDP_VALUE(int, n) IDP_GRID) {
  int i = IDP_GROUP * IDP_WIDTH + IDP_LANE;
  if (i < n)
    out[i] = 0;
}
IDP_KERNEL void idp_fraction_kernel(IDP_GLOBAL real *fraction, IDP_GLOBAL const real *coeff,
                                    IDP_GLOBAL const real *h, IDP_GLOBAL const real *wave,
                                    IDP_VALUE(real, c), IDP_VALUE(int, n) IDP_GRID) {
  int i = IDP_GROUP * IDP_WIDTH + IDP_LANE;
  if (i < n) {
    real low = c * h[i] * wave[i];
    fraction[i] = low > 0 ? rmin(1.0f, rmax(0.0f, coeff[i] / low)) : 0;
  }
}
#endif
IDP_INLINE void load_state(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                           IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                           IDP_GLOBAL const real *energy, int i, IDP_PRIVATE real s[5]) {
  s[0] = rho[i];
  s[1] = mx[i];
  s[2] = my[i];
  s[3] = mz[i];
  s[4] = energy[i];
}

IDP_INLINE real internal_energy(IDP_PRIVATE const real s[5]) {
  return s[4] - 0.5f * (s[1] * s[1] + s[2] * s[2] + s[3] * s[3]) / s[0];
}

IDP_INLINE bool state_is_admissible(IDP_PRIVATE const real s[5], real floor) {
  for (int c = 0; c < 5; ++c)
    if (!isfinite(s[c]))
      return false;
  if (s[0] <= 0.0f)
    return false;
  const real margin =
      s[0] * (s[4] - floor) - 0.5f * (s[1] * s[1] + s[2] * s[2] + s[3] * s[3]);
  return isfinite(margin) && margin >= 0.0f;
}

IDP_INLINE real specific_entropy(IDP_PRIVATE const real s[5], real gamma) {
  const real huge = IDP_HUGE;
  for (int c = 0; c < 5; ++c)
    if (!isfinite(s[c]))
      return -huge;
  if (s[0] <= 0.0f)
    return -huge;
  const real pressure = (gamma - 1.0f) * internal_energy(s);
  if (!isfinite(pressure) || pressure <= 0.0f)
    return -huge;
  return log(pressure) - gamma * log(s[0]);
}

IDP_INLINE real entropy_tolerance(IDP_PRIVATE const real s[5], real lower) {
  const real eps = IDP_EPS;
  const real tiny = IDP_TINY;
  real conditioning = 1.0f;
  bool finite = s[0] > 0.0f;
  for (int c = 0; c < 5; ++c)
    finite = finite && isfinite(s[c]);
  if (finite) {
    const real kinetic = 0.5f * (s[1] * s[1] + s[2] * s[2] + s[3] * s[3]) / s[0];
    const real internal = s[4] - kinetic;
    if (isfinite(internal) && internal > tiny)
      conditioning = (fabs(s[4]) + fabs(kinetic)) / internal;
  }
  return rmin(sqrt(eps),
              64.0f * eps * rmax(1.0f, rmax(fabs(lower), conditioning)));
}

IDP_INLINE bool entropy_is_admissible(IDP_PRIVATE const real s[5], real gamma, real lower,
                                      real floor) {
  if (!state_is_admissible(s, floor))
    return false;
  const real entropy = specific_entropy(s, gamma);
  return isfinite(entropy) && entropy >= lower - entropy_tolerance(s, lower);
}

IDP_INLINE void fluxes(IDP_PRIVATE const real s[5], real gamma, IDP_PRIVATE real fx[5],
                       IDP_PRIVATE real fy[5], IDP_PRIVATE real fz[5]) {
  const real inv_rho = 1.0f / s[0];
  const real u = s[1] * inv_rho;
  const real v = s[2] * inv_rho;
  const real w = s[3] * inv_rho;
  const real p = (gamma - 1.0f) * internal_energy(s);
  fx[0] = s[1];
  fy[0] = s[2];
  fz[0] = s[3];
  fx[1] = s[1] * u + p;
  fy[1] = s[1] * v;
  fz[1] = s[1] * w;
  fx[2] = s[2] * u;
  fy[2] = s[2] * v + p;
  fz[2] = s[2] * w;
  fx[3] = s[3] * u;
  fy[3] = s[3] * v;
  fz[3] = s[3] * w + p;
  fx[4] = (s[4] + p) * u;
  fy[4] = (s[4] + p) * v;
  fz[4] = (s[4] + p) * w;
}

IDP_INLINE real ordered_wave_speed(IDP_PRIVATE const real left[5],
                                   IDP_PRIVATE const real right[5],
                                   IDP_PRIVATE const real normal[3], real gamma) {
  const real pl = (gamma - 1.0f) * internal_energy(left);
  const real pr = (gamma - 1.0f) * internal_energy(right);
  const real al = sqrt(gamma * pl / left[0]);
  const real ar = sqrt(gamma * pr / right[0]);
  const real ul = (left[1] * normal[0] + left[2] * normal[1] + left[3] * normal[2]) / left[0];
  const real ur =
      (right[1] * normal[0] + right[2] * normal[1] + right[3] * normal[2]) / right[0];
  real pmin, pmax, dmin, amin, amax;
  if (pl <= pr) {
    pmin = pl;
    pmax = pr;
    dmin = left[0];
    amin = al;
    amax = ar;
  } else {
    pmin = pr;
    pmax = pl;
    dmin = right[0];
    amin = ar;
    amax = al;
  }
  const real exponent = (gamma - 1.0f) / (2.0f * gamma);
  const real ratio = pow(pmin / pmax, exponent);
  const real phi_min =
      2.0f * amax * (ratio - 1.0f) / (gamma - 1.0f) + ur - ul;
  if (phi_min >= 0.0f)
    return rmax(rmax(-(ul - al), 0.0f), rmax(ur + ar, 0.0f));
  const real ca = 2.0f / ((gamma + 1.0f) * dmin);
  const real cb = pmin * (gamma - 1.0f) / (gamma + 1.0f);
  const real phi_max = (pmax - pmin) * sqrt(ca / (pmax + cb)) + ur - ul;
  const real numerator = amin + amax - 0.5f * (gamma - 1.0f) * (ur - ul);
  const real p_two = pmin * pow(numerator / (amin + amax * ratio),
                                2.0f * gamma / (gamma - 1.0f));
  const real p_upper = phi_max < 0.0f ? p_two : rmin(pmax, p_two);
  const real sl =
      ul - al * sqrt(1.0f + rmax((p_upper - pl) / pl, 0.0f) *
                                       (gamma + 1.0f) / (2.0f * gamma));
  const real sr =
      ur + ar * sqrt(1.0f + rmax((p_upper - pr) / pr, 0.0f) *
                                       (gamma + 1.0f) / (2.0f * gamma));
  return rmax(rmax(-sl, 0.0f), rmax(sr, 0.0f));
}

IDP_INLINE real maximum_wave_speed(IDP_PRIVATE const real left[5],
                                   IDP_PRIVATE const real right[5],
                                   IDP_PRIVATE const real normal[3], real gamma) {
  real reverse[3] = {-normal[0], -normal[1], -normal[2]};
  return rmax(ordered_wave_speed(left, right, normal, gamma),
              ordered_wave_speed(right, left, reverse, gamma));
}

IDP_INLINE real energy_margin(IDP_PRIVATE const real s[5], real floor) {
  return s[4] - floor - 0.5f * (s[1] * s[1] + s[2] * s[2] + s[3] * s[3]) / s[0];
}

IDP_INLINE real energy_margin_derivative(IDP_PRIVATE const real s[5],
                                         IDP_PRIVATE const real q[5]) {
  const real m2 = s[1] * s[1] + s[2] * s[2] + s[3] * s[3];
  return q[4] - (s[1] * q[1] + s[2] * q[2] + s[3] * q[3]) / s[0] +
         0.5f * m2 * q[0] / (s[0] * s[0]);
}

IDP_INLINE real entropy_margin(IDP_PRIVATE const real s[5], real lower, real gamma) {
  return (gamma - 1.0f) * energy_margin(s, 0.0f) - exp(lower) * pow(s[0], gamma);
}

IDP_INLINE real entropy_margin_derivative(IDP_PRIVATE const real s[5],
                                          IDP_PRIVATE const real q[5], real lower, real gamma) {
  return (gamma - 1.0f) * energy_margin_derivative(s, q) -
         exp(lower) * gamma * pow(s[0], gamma - 1.0f) * q[0];
}

IDP_INLINE void limit_energy(IDP_PRIVATE const real base[5],
                             IDP_PRIVATE const real correction[5], real floor,
                             IDP_PRIVATE real *limit) {
  const real eps = IDP_EPS;
  real left = 0.0f, right = (*limit), trial[5];
  real vl = energy_margin(base, floor);
  for (int c = 0; c < 5; ++c)
    trial[c] = base[c] + right * correction[c];
  real vr = energy_margin(trial, floor);
  const real tol = 256.0f * eps * rmax(1.0f, rmax(fabs(vl), fabs(vr)));
  if (vl <= 0.0f) {
    (*limit) = 0.0f;
    return;
  }
  for (int it = 0; it < 16; ++it) {
    if (right - left <= 256.0f * eps * rmax(1.0f, right) || vl <= vr)
      break;
    const real slope = (vr - vl) / (right - left);
    real next = left - vl / slope;
    if (!isfinite(next) || next <= left || next >= right)
      break;
    for (int c = 0; c < 5; ++c)
      trial[c] = base[c] + next * correction[c];
    real vn = energy_margin(trial, floor);
    if (!isfinite(vn) || vn < -tol)
      break;
    left = next;
    vl = vn;
    if (fabs(vn) <= tol)
      break;
    for (int c = 0; c < 5; ++c)
      trial[c] = base[c] + right * correction[c];
    const real ds = energy_margin_derivative(trial, correction);
    if (!isfinite(ds) || ds >= 0.0f)
      break;
    next = right - vr / ds;
    if (!isfinite(next) || next <= left || next >= right)
      break;
    for (int c = 0; c < 5; ++c)
      trial[c] = base[c] + next * correction[c];
    vn = energy_margin(trial, floor);
    if (!isfinite(vn) || vn > tol)
      break;
    if (vn >= 0.0f) {
      left = next;
      break;
    }
    right = next;
    vr = vn;
  }
  (*limit) = left;
}

IDP_INLINE void limit_entropy(IDP_PRIVATE const real base[5],
                              IDP_PRIVATE const real correction[5], real lower, real gamma,
                              IDP_PRIVATE real *limit) {
  const real eps = IDP_EPS;
  real left = 0.0f, right = (*limit), trial[5];
  real vl = entropy_margin(base, lower, gamma);
  for (int c = 0; c < 5; ++c)
    trial[c] = base[c] + right * correction[c];
  real vr = entropy_margin(trial, lower, gamma);
  const real tol = 256.0f * eps * rmax(1.0f, rmax(fabs(vl), fabs(vr)));
  if (vl <= 0.0f) {
    (*limit) = 0.0f;
    return;
  }
  for (int it = 0; it < 16; ++it) {
    if (right - left <= 256.0f * eps * rmax(1.0f, right) || vl <= vr)
      break;
    const real slope = (vr - vl) / (right - left);
    real next = left - vl / slope;
    if (!isfinite(next) || next <= left || next >= right)
      break;
    for (int c = 0; c < 5; ++c)
      trial[c] = base[c] + next * correction[c];
    real vn = entropy_margin(trial, lower, gamma);
    if (!isfinite(vn) || vn < -tol)
      break;
    left = next;
    vl = vn;
    if (fabs(vn) <= tol)
      break;
    for (int c = 0; c < 5; ++c)
      trial[c] = base[c] + right * correction[c];
    const real ds = entropy_margin_derivative(trial, correction, lower, gamma);
    if (!isfinite(ds) || ds >= 0.0f)
      break;
    next = right - vr / ds;
    if (!isfinite(next) || next <= left || next >= right)
      break;
    for (int c = 0; c < 5; ++c)
      trial[c] = base[c] + next * correction[c];
    vn = entropy_margin(trial, lower, gamma);
    if (!isfinite(vn) || vn > tol)
      break;
    if (vn >= 0.0f) {
      left = next;
      break;
    }
    right = next;
    vr = vn;
  }
  (*limit) = left;
}

IDP_INLINE void entropy_bisection(IDP_PRIVATE const real base[5],
                                  IDP_PRIVATE const real correction[5], real lower_entropy,
                                  real gamma, real floor, IDP_PRIVATE real *limit) {
  real trial[5];
  for (int c = 0; c < 5; ++c)
    trial[c] = base[c] + (*limit) * correction[c];
  if (entropy_is_admissible(trial, gamma, lower_entropy, floor))
    return;
  if (!entropy_is_admissible(base, gamma, lower_entropy, floor)) {
    (*limit) = 0.0f;
    return;
  }
  real left = 0.0f, right = (*limit);
  for (int it = 0; it < 64; ++it) {
    const real mid = left + 0.5f * (right - left);
    if (mid <= left || mid >= right)
      break;
    for (int c = 0; c < 5; ++c)
      trial[c] = base[c] + mid * correction[c];
    if (entropy_is_admissible(trial, gamma, lower_entropy, floor))
      left = mid;
    else
      right = mid;
  }
  (*limit) = left;
}

IDP_INLINE real limit_endpoint(IDP_PRIVATE const real base[5],
                               IDP_PRIVATE const real correction[5], real dl, real du, real sl,
                               real gamma, real floor, bool enforce_energy,
                               bool enforce_entropy, bool check_base,
                               IDP_PRIVATE bool *density_limited,
                               IDP_PRIVATE bool *energy_limited,
                               IDP_PRIVATE bool *entropy_limited) {
  const real eps = IDP_EPS;
  real limit = 1.0f, trial[5];
  (*density_limited) = false;
  (*energy_limited) = false;
  (*entropy_limited) = false;
  if (check_base) {
    bool finite = true;
    for (int c = 0; c < 5; ++c)
      finite = finite && isfinite(base[c]);
    if (!finite) {
      (*energy_limited) = enforce_energy;
      (*entropy_limited) = enforce_entropy;
      return 0.0f;
    }
    if (base[0] <= 0.0f || base[0] < dl || base[0] > du) {
      (*density_limited) = true;
      return 0.0f;
    }
    if (enforce_energy && !state_is_admissible(base, floor)) {
      (*energy_limited) = true;
      return 0.0f;
    }
    if (enforce_entropy && !entropy_is_admissible(base, gamma, sl, floor)) {
      (*entropy_limited) = true;
      return 0.0f;
    }
  }
  if (correction[0] < 0.0f)
    limit = rmin(limit, (base[0] - dl) / (-correction[0]));
  else if (correction[0] > 0.0f)
    limit = rmin(limit, (du - base[0]) / correction[0]);
  if (correction[0] != 0.0f) {
    limit = rmax(0.0f, limit);
    (*density_limited) = limit < 1.0f;
    if ((*density_limited) && limit > 0.0f)
      limit *= 1.0f - 32.0f * eps;
  }
  for (int c = 0; c < 5; ++c)
    trial[c] = base[c] + limit * correction[c];
  (*energy_limited) = enforce_energy && !state_is_admissible(trial, floor);
  (*entropy_limited) = enforce_entropy && !entropy_is_admissible(trial, gamma, sl, floor);
  if (!(*energy_limited) && !(*entropy_limited))
    return limit;
  if ((*energy_limited))
    limit_energy(base, correction, floor, &limit);
  for (int c = 0; c < 5; ++c)
    trial[c] = base[c] + limit * correction[c];
  if (enforce_entropy && !entropy_is_admissible(trial, gamma, sl, floor)) {
    limit_entropy(base, correction, sl, gamma, &limit);
    (*entropy_limited) = true;
  }
  if (limit > 0.0f)
    limit *= 1.0f - 32.0f * eps;
  if (enforce_entropy && limit > 0.0f)
    entropy_bisection(base, correction, sl, gamma, floor, &limit);
  return limit;
}

IDP_INLINE real floor_timestep(IDP_PRIVATE const real state[5],
                               IDP_PRIVATE const real residual[5], real floor, real upper) {
  const real eps = IDP_EPS;
  real limit = rmax(0.0f, upper), trial[5];
  if (limit <= 0.0f)
    return limit;
  if (residual[0] > 0.0f) {
    const real density_limit = state[0] / residual[0];
    if (density_limit <= limit)
      limit = rmax(0.0f, density_limit * (1.0f - 64.0f * eps));
  }
  if (limit <= 0.0f)
    return limit;
  for (int c = 0; c < 5; ++c)
    trial[c] = state[c] - limit * residual[c];
  if (state_is_admissible(trial, floor))
    return limit;
  if (!state_is_admissible(state, floor))
    return 0.0f;
  // Leave room for roundoff when the caller forms the candidate again.
  const real guarded_floor =
      floor + 16.0f * eps * rmax(fabs(state[4]), limit * fabs(residual[4]));
  real left = 0.0f, right = limit;
  for (int it = 0; it < 64; ++it) {
    const real mid = 0.5f * (left + right);
    for (int c = 0; c < 5; ++c)
      trial[c] = state[c] - mid * residual[c];
    if (state_is_admissible(trial, guarded_floor))
      left = mid;
    else
      right = mid;
  }
  return left;
}

IDP_KERNEL void idp_primitives_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                                      IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                                      IDP_GLOBAL const real *energy, IDP_GLOBAL real *u,
                                      IDP_GLOBAL real *v, IDP_GLOBAL real *w,
                                      IDP_GLOBAL real *p, IDP_GLOBAL real *sound,
                                      IDP_GLOBAL real *internal, IDP_GLOBAL real *status,
                                      IDP_VALUE(real, gamma), IDP_VALUE(real, floor),
                                      IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    real s[5];
    load_state(rho, mx, my, mz, energy, i, s);
    status[i] = 0.0f;
    u[i] = v[i] = w[i] = p[i] = sound[i] = internal[i] = 0.0f;
    if (!isfinite(gamma) || gamma <= 1.0f || !state_is_admissible(s, floor)) {
      status[i] = 1.0f;
      continue;
    }
    u[i] = s[1] / s[0];
    v[i] = s[2] / s[0];
    w[i] = s[3] / s[0];
    internal[i] = internal_energy(s);
    p[i] = (gamma - 1.0f) * internal[i];
    sound[i] = sqrt(gamma * p[i] / s[0]);
    if (!isfinite(sound[i]))
      status[i] = 1.0f;
  }
}

IDP_KERNEL void idp_graph_viscosity_kernel(
    IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx, IDP_GLOBAL const real *my,
    IDP_GLOBAL const real *mz, IDP_GLOBAL const real *energy, IDP_GLOBAL const real *user_speed,
    IDP_VALUE(int, has_user), IDP_GLOBAL const int *left, IDP_GLOBAL const int *right,
    IDP_GLOBAL const real *coef, IDP_GLOBAL real *edge_visc, IDP_GLOBAL real *vis_sum,
    IDP_GLOBAL real *edge_speed, IDP_VALUE(real, gamma), IDP_VALUE(int, nedge) IDP_GRID) {
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nedge; e += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[e], b = right[e];
    const real cx = coef[3 * e], cy = coef[3 * e + 1], cz = coef[3 * e + 2];
    const real norm = sqrt(cx * cx + cy * cy + cz * cz);
    real wave;
    if (has_user)
      wave = rmax(user_speed[a], user_speed[b]);
    else {
      IDP_PRIVATE real l[5], r[5], normal[3] = {cx / norm, cy / norm, cz / norm};
      load_state(rho, mx, my, mz, energy, a, l);
      load_state(rho, mx, my, mz, energy, b, r);
      wave = maximum_wave_speed(l, r, normal, gamma);
    }
    edge_speed[e] = wave;
    edge_visc[e] = norm * wave;
    idp_add(vis_sum + a, edge_visc[e]);
    idp_add(vis_sum + b, edge_visc[e]);
  }
}

IDP_KERNEL void idp_graph_rate_kernel(IDP_GLOBAL const real *vis, IDP_GLOBAL const real *mass,
                                      IDP_GLOBAL real *rate, IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS)
    rate[i] = 2.0f * vis[i] / mass[i];
}

IDP_KERNEL void idp_nodal_wave_kernel(IDP_GLOBAL const real *u, IDP_GLOBAL const real *v,
                                      IDP_GLOBAL const real *w, IDP_GLOBAL const real *sound,
                                      IDP_GLOBAL real *wave, IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS)
    wave[i] = sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]) + sound[i];
}

IDP_KERNEL void idp_low_init_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                                    IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                                    IDP_GLOBAL const real *energy, IDP_GLOBAL const real *diag,
                                    IDP_GLOBAL real *r0, IDP_GLOBAL real *r1,
                                    IDP_GLOBAL real *r2, IDP_GLOBAL real *r3,
                                    IDP_GLOBAL real *r4, IDP_VALUE(real, gamma),
                                    IDP_VALUE(int, periodic), IDP_VALUE(int, n) IDP_GRID) {
  IDP_GLOBAL real *out[5] = {r0, r1, r2, r3, r4};
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    if (periodic) {
      for (int c = 0; c < 5; ++c)
        out[c][i] = 0.0f;
      continue;
    }
    real s[5], fx[5], fy[5], fz[5];
    load_state(rho, mx, my, mz, energy, i, s);
    fluxes(s, gamma, fx, fy, fz);
    for (int c = 0; c < 5; ++c)
      out[c][i] = diag[3 * i] * fx[c] + diag[3 * i + 1] * fy[c] + diag[3 * i + 2] * fz[c];
  }
}

IDP_KERNEL void idp_low_edge_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                                    IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                                    IDP_GLOBAL const real *energy, IDP_GLOBAL const int *left,
                                    IDP_GLOBAL const int *right, IDP_GLOBAL const real *coef,
                                    IDP_GLOBAL const real *visc, IDP_GLOBAL real *r0,
                                    IDP_GLOBAL real *r1, IDP_GLOBAL real *r2,
                                    IDP_GLOBAL real *r3, IDP_GLOBAL real *r4,
                                    IDP_VALUE(real, gamma), IDP_VALUE(int, periodic),
                                    IDP_VALUE(int, nedge) IDP_GRID) {
  IDP_GLOBAL real *out[5] = {r0, r1, r2, r3, r4};
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nedge; e += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[e], b = right[e];
    real sa[5], sb[5], fax[5], fay[5], faz[5];
    real fbx[5], fby[5], fbz[5];
    load_state(rho, mx, my, mz, energy, a, sa);
    load_state(rho, mx, my, mz, energy, b, sb);
    fluxes(sa, gamma, fax, fay, faz);
    fluxes(sb, gamma, fbx, fby, fbz);
    const real cx = coef[3 * e], cy = coef[3 * e + 1], cz = coef[3 * e + 2];
    for (int c = 0; c < 5; ++c) {
      const real fa = cx * fax[c] + cy * fay[c] + cz * faz[c];
      const real fb = cx * fbx[c] + cy * fby[c] + cz * fbz[c];
      const real d = visc[e] * (sb[c] - sa[c]);
      if (periodic) {
        idp_add(out[c] + a, fb - fa - d);
        idp_add(out[c] + b, fb - fa + d);
      } else {
        idp_add(out[c] + a, fb - d);
        idp_add(out[c] + b, -fa + d);
      }
    }
  }
}

IDP_KERNEL void idp_scale_residual_kernel(IDP_GLOBAL real *r0, IDP_GLOBAL real *r1,
                                          IDP_GLOBAL real *r2, IDP_GLOBAL real *r3,
                                          IDP_GLOBAL real *r4, IDP_GLOBAL const real *binv,
                                          IDP_VALUE(int, n) IDP_GRID) {
  IDP_GLOBAL real *out[5] = {r0, r1, r2, r3, r4};
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS)
    for (int c = 0; c < 5; ++c)
      out[c][i] *= binv[i];
}

IDP_KERNEL void idp_correction_residual_init_kernel(
    IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx, IDP_GLOBAL const real *my,
    IDP_GLOBAL const real *mz, IDP_GLOBAL const real *energy, IDP_GLOBAL const real *diag,
    IDP_GLOBAL const real *drdx, IDP_GLOBAL const real *drdy, IDP_GLOBAL const real *drdz,
    IDP_GLOBAL const real *dsdx, IDP_GLOBAL const real *dsdy, IDP_GLOBAL const real *dsdz,
    IDP_GLOBAL const real *dtdx, IDP_GLOBAL const real *dtdy, IDP_GLOBAL const real *dtdz,
    IDP_GLOBAL const real *dx, IDP_GLOBAL const real *dy, IDP_GLOBAL const real *dz,
    IDP_GLOBAL const real *wx, IDP_GLOBAL const real *wy, IDP_GLOBAL const real *wz,
    IDP_GLOBAL real *r0, IDP_GLOBAL real *r1, IDP_GLOBAL real *r2, IDP_GLOBAL real *r3,
    IDP_GLOBAL real *r4, IDP_VALUE(real, gamma), IDP_VALUE(int, lx), IDP_VALUE(int, ly),
    IDP_VALUE(int, lz), IDP_VALUE(int, n) IDP_GRID) {
  IDP_GLOBAL real *out[5] = {r0, r1, r2, r3, r4};
  for (int node = IDP_GROUP * IDP_WIDTH + IDP_LANE; node < n; node += IDP_WIDTH * IDP_GROUPS) {
    int q = node;
    const int i = q % lx;
    q /= lx;
    const int j = q % ly;
    q /= ly;
    const int k = q % lz;
    const int e = q / lz;
    real s[5], fx[5], fy[5], fz[5];
    load_state(rho, mx, my, mz, energy, node, s);
    fluxes(s, gamma, fx, fy, fz);
    real high[5] = {0, 0, 0, 0, 0};
    const real weight = wx[i] * wy[j] * wz[k];
    for (int l = 0; l < lx; ++l) {
      const int nl = ((e * lz + k) * ly + j) * lx + l;
      real sl[5], x[5], y[5], z[5];
      load_state(rho, mx, my, mz, energy, nl, sl);
      fluxes(sl, gamma, x, y, z);
      const real d = weight * dx[i + l * lx];
      for (int c = 0; c < 5; ++c)
        high[c] += d * (drdx[nl] * x[c] + drdy[nl] * y[c] + drdz[nl] * z[c]);
    }
    for (int l = 0; l < ly; ++l) {
      const int nl = ((e * lz + k) * ly + l) * lx + i;
      real sl[5], x[5], y[5], z[5];
      load_state(rho, mx, my, mz, energy, nl, sl);
      fluxes(sl, gamma, x, y, z);
      const real d = weight * dy[j + l * ly];
      for (int c = 0; c < 5; ++c)
        high[c] += d * (dsdx[nl] * x[c] + dsdy[nl] * y[c] + dsdz[nl] * z[c]);
    }
    if (lz > 1)
      for (int l = 0; l < lz; ++l) {
        const int nl = ((e * lz + l) * ly + j) * lx + i;
        real sl[5], x[5], y[5], z[5];
        load_state(rho, mx, my, mz, energy, nl, sl);
        fluxes(sl, gamma, x, y, z);
        const real d = weight * dz[k + l * lz];
        for (int c = 0; c < 5; ++c)
          high[c] += d * (dtdx[nl] * x[c] + dtdy[nl] * y[c] + dtdz[nl] * z[c]);
      }
    for (int c = 0; c < 5; ++c)
      out[c][node] = diag[3 * node] * fx[c] + diag[3 * node + 1] * fy[c] +
                     diag[3 * node + 2] * fz[c] - high[c];
  }
}

IDP_KERNEL void idp_correction_residual_edge_kernel(
    IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx, IDP_GLOBAL const real *my,
    IDP_GLOBAL const real *mz, IDP_GLOBAL const real *energy, IDP_GLOBAL const int *left,
    IDP_GLOBAL const int *right, IDP_GLOBAL const real *coef, IDP_GLOBAL real *r0,
    IDP_GLOBAL real *r1, IDP_GLOBAL real *r2, IDP_GLOBAL real *r3, IDP_GLOBAL real *r4,
    IDP_VALUE(real, gamma), IDP_VALUE(int, nedge) IDP_GRID) {
  IDP_GLOBAL real *out[5] = {r0, r1, r2, r3, r4};
  for (int edge = IDP_GROUP * IDP_WIDTH + IDP_LANE; edge < nedge;
       edge += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[edge], b = right[edge];
    real sa[5], sb[5], fax[5], fay[5], faz[5];
    real fbx[5], fby[5], fbz[5];
    load_state(rho, mx, my, mz, energy, a, sa);
    load_state(rho, mx, my, mz, energy, b, sb);
    fluxes(sa, gamma, fax, fay, faz);
    fluxes(sb, gamma, fbx, fby, fbz);
    const real cx = coef[3 * edge];
    const real cy = coef[3 * edge + 1], cz = coef[3 * edge + 2];
    for (int c = 0; c < 5; ++c) {
      idp_add(out[c] + a, cx * fbx[c] + cy * fby[c] + cz * fbz[c]);
      idp_add(out[c] + b, -cx * fax[c] - cy * fay[c] - cz * faz[c]);
    }
  }
}

IDP_KERNEL void idp_reconstruct_general_kernel(
    IDP_GLOBAL const real *r0, IDP_GLOBAL const real *r1, IDP_GLOBAL const real *r2,
    IDP_GLOBAL const real *r3, IDP_GLOBAL const real *r4, IDP_GLOBAL const real *wx,
    IDP_GLOBAL const real *wy, IDP_GLOBAL const real *wz, IDP_GLOBAL real *correction,
    IDP_GLOBAL real *error, IDP_VALUE(int, lx), IDP_VALUE(int, ly), IDP_VALUE(int, lz),
    IDP_VALUE(int, edges_per_element), IDP_VALUE(int, nelv) IDP_GRID) {
  IDP_GLOBAL const real *res[5] = {r0, r1, r2, r3, r4};
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nelv; e += IDP_WIDTH * IDP_GROUPS) {
    real sx = 0.0f, sy = 0.0f, sz = 0.0f;
    for (int i = 0; i < lx; ++i)
      sx += wx[i];
    for (int j = 0; j < ly; ++j)
      sy += wy[j];
    for (int k = 0; k < lz; ++k)
      sz += wz[k];
    real element_sum[5] = {0, 0, 0, 0, 0}, scale[5] = {0, 0, 0, 0, 0};
    for (int k = 0; k < lz; ++k)
      for (int j = 0; j < ly; ++j)
        for (int i = 0; i < lx; ++i) {
          const int node = ((e * lz + k) * ly + j) * lx + i;
          for (int c = 0; c < 5; ++c) {
            element_sum[c] += res[c][node];
            scale[c] += fabs(res[c][node]);
          }
        }
    for (int c = 0; c < 5; ++c)
      idp_maximum(error + c, fabs(element_sum[c]) / rmax(1.0f, scale[c]));
    const int base = e * edges_per_element;
    int edge = base;
    for (int k = 0; k < lz; ++k)
      for (int j = 0; j < ly; ++j) {
        real line[5] = {0, 0, 0, 0, 0}, cum[5] = {0, 0, 0, 0, 0};
        for (int i = 0; i < lx; ++i) {
          const int n = ((e * lz + k) * ly + j) * lx + i;
          for (int c = 0; c < 5; ++c)
            line[c] += res[c][n];
        }
        for (int i = 0; i < lx - 1; ++i) {
          const int n = ((e * lz + k) * ly + j) * lx + i;
          for (int c = 0; c < 5; ++c) {
            cum[c] += res[c][n] - wx[i] / sx * line[c];
            correction[5 * edge + c] = cum[c];
          }
          ++edge;
        }
      }
    const int ybase = base + (lx - 1) * ly * lz;
    edge = ybase;
    for (int k = 0; k < lz; ++k) {
      real plane[5] = {0, 0, 0, 0, 0};
      for (int j = 0; j < ly; ++j)
        for (int i = 0; i < lx; ++i) {
          const int n = ((e * lz + k) * ly + j) * lx + i;
          for (int c = 0; c < 5; ++c)
            plane[c] += res[c][n];
        }
      for (int i = 0; i < lx; ++i) {
        real cum[5] = {0, 0, 0, 0, 0};
        for (int j = 0; j < ly - 1; ++j) {
          real line[5] = {0, 0, 0, 0, 0};
          for (int ii = 0; ii < lx; ++ii) {
            const int n = ((e * lz + k) * ly + j) * lx + ii;
            for (int c = 0; c < 5; ++c)
              line[c] += res[c][n];
          }
          for (int c = 0; c < 5; ++c) {
            cum[c] += wx[i] / sx * (line[c] - wy[j] / sy * plane[c]);
            correction[5 * edge + c] = cum[c];
          }
          ++edge;
        }
      }
    }
    edge = ybase + lx * (ly - 1) * lz;
    if (lz > 1)
      for (int j = 0; j < ly; ++j)
        for (int i = 0; i < lx; ++i) {
          real cum[5] = {0, 0, 0, 0, 0};
          for (int k = 0; k < lz - 1; ++k) {
            real plane[5] = {0, 0, 0, 0, 0};
            for (int jj = 0; jj < ly; ++jj)
              for (int ii = 0; ii < lx; ++ii) {
                const int n = ((e * lz + k) * ly + jj) * lx + ii;
                for (int c = 0; c < 5; ++c)
                  plane[c] += res[c][n];
              }
            for (int c = 0; c < 5; ++c) {
              cum[c] += wx[i] / sx * wy[j] / sy * (plane[c] - wz[k] / sz * element_sum[c]);
              correction[5 * edge + c] = cum[c];
            }
            ++edge;
          }
        }
  }
}

IDP_INLINE void line_node_residual(int direction, int pos, int i, int j, int k, int e, int lx,
                                   int ly, int lz, int edge_first, IDP_GLOBAL const real *rho,
                                   IDP_GLOBAL const real *mx, IDP_GLOBAL const real *my,
                                   IDP_GLOBAL const real *mz, IDP_GLOBAL const real *energy,
                                   IDP_GLOBAL const real *coef, IDP_GLOBAL const real *B,
                                   IDP_GLOBAL const real *jacinv, IDP_GLOBAL const real *mr0,
                                   IDP_GLOBAL const real *mr1, IDP_GLOBAL const real *mr2,
                                   IDP_GLOBAL const real *D, real gamma,
                                   IDP_PRIVATE real low[5], IDP_PRIVATE real high[5]) {
  const int length = direction == 0 ? lx : (direction == 1 ? ly : lz);
  int ci = i, cj = j, ck = k;
  if (direction == 0)
    ci = pos;
  else if (direction == 1)
    cj = pos;
  else
    ck = pos;
  const int node = ((e * lz + ck) * ly + cj) * lx + ci;
  for (int c = 0; c < 5; ++c) {
    low[c] = 0.0f;
    high[c] = 0.0f;
  }
  if (pos < length - 1) {
    int ni = ci, nj = cj, nk = ck;
    if (direction == 0)
      ++ni;
    else if (direction == 1)
      ++nj;
    else
      ++nk;
    int nn = ((e * lz + nk) * ly + nj) * lx + ni;
    real s[5], x[5], y[5], z[5];
    load_state(rho, mx, my, mz, energy, nn, s);
    fluxes(s, gamma, x, y, z);
    int ed = edge_first + pos;
    for (int c = 0; c < 5; ++c)
      low[c] += coef[3 * ed] * x[c] + coef[3 * ed + 1] * y[c] + coef[3 * ed + 2] * z[c];
  }
  if (pos > 0) {
    int ni = ci, nj = cj, nk = ck;
    if (direction == 0)
      --ni;
    else if (direction == 1)
      --nj;
    else
      --nk;
    int nn = ((e * lz + nk) * ly + nj) * lx + ni;
    real s[5], x[5], y[5], z[5];
    load_state(rho, mx, my, mz, energy, nn, s);
    fluxes(s, gamma, x, y, z);
    int ed = edge_first + pos - 1;
    for (int c = 0; c < 5; ++c)
      low[c] -= coef[3 * ed] * x[c] + coef[3 * ed + 1] * y[c] + coef[3 * ed + 2] * z[c];
  }
  if (pos == 0 || pos == length - 1) {
    real s[5], x[5], y[5], z[5];
    load_state(rho, mx, my, mz, energy, node, s);
    fluxes(s, gamma, x, y, z);
    int ed = pos == 0 ? edge_first : edge_first + length - 2;
    const real sign = pos == 0 ? -1.0f : 1.0f;
    for (int c = 0; c < 5; ++c)
      low[c] +=
          sign * (coef[3 * ed] * x[c] + coef[3 * ed + 1] * y[c] + coef[3 * ed + 2] * z[c]);
  }
  for (int l = 0; l < length; ++l) {
    int ni = i, nj = j, nk = k;
    if (direction == 0)
      ni = l;
    else if (direction == 1)
      nj = l;
    else
      nk = l;
    int nn = ((e * lz + nk) * ly + nj) * lx + ni;
    real s[5], x[5], y[5], z[5];
    load_state(rho, mx, my, mz, energy, nn, s);
    fluxes(s, gamma, x, y, z);
    const real d = D[pos + l * length];
    for (int c = 0; c < 5; ++c)
      high[c] += d * (mr0[node] * x[c] + mr1[node] * y[c] + mr2[node] * z[c]);
  }
  for (int c = 0; c < 5; ++c)
    high[c] *= B[node] * jacinv[node];
}

IDP_KERNEL void idp_reconstruct_direction_kernel(
    IDP_VALUE(int, direction), IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
    IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz, IDP_GLOBAL const real *energy,
    IDP_GLOBAL const real *coef, IDP_GLOBAL const real *B, IDP_GLOBAL const real *jacinv,
    IDP_GLOBAL const real *mr0, IDP_GLOBAL const real *mr1, IDP_GLOBAL const real *mr2,
    IDP_GLOBAL const real *D, IDP_GLOBAL real *correction, IDP_GLOBAL real *error,
    IDP_VALUE(real, gamma), IDP_VALUE(int, lx), IDP_VALUE(int, ly), IDP_VALUE(int, lz),
    IDP_VALUE(int, edges_per_element), IDP_VALUE(int, line_count) IDP_GRID) {
  for (int line = IDP_GROUP * IDP_WIDTH + IDP_LANE; line < line_count;
       line += IDP_WIDTH * IDP_GROUPS) {
    int q = line, i = 0, j = 0, k = 0, e = 0, first = 0, length = 0;
    const int xcount = (lx - 1) * ly * lz, ycount = lx * (ly - 1) * lz;
    if (direction == 0) {
      j = q % ly;
      q /= ly;
      k = q % lz;
      e = q / lz;
      length = lx;
      first = e * edges_per_element + (k * ly + j) * (lx - 1);
    } else if (direction == 1) {
      i = q % lx;
      q /= lx;
      k = q % lz;
      e = q / lz;
      length = ly;
      first = e * edges_per_element + xcount + (k * lx + i) * (ly - 1);
    } else {
      i = q % lx;
      q /= lx;
      j = q % ly;
      e = q / ly;
      length = lz;
      first = e * edges_per_element + xcount + ycount + (j * lx + i) * (lz - 1);
    }
    real cumulative[5] = {0, 0, 0, 0, 0}, sum[5] = {0, 0, 0, 0, 0}, scale[5] = {0, 0, 0, 0, 0};
    for (int pos = 0; pos < length; ++pos) {
      real low[5], high[5];
      line_node_residual(direction, pos, i, j, k, e, lx, ly, lz, first, rho, mx, my, mz, energy,
                         coef, B, jacinv, mr0, mr1, mr2, D, gamma, low, high);
      for (int c = 0; c < 5; ++c) {
        const real value = low[c] - high[c];
        sum[c] += value;
        scale[c] += fabs(low[c]) + fabs(high[c]);
        if (pos < length - 1) {
          cumulative[c] += value;
          correction[5 * (first + pos) + c] = cumulative[c];
        }
      }
    }
    for (int c = 0; c < 5; ++c)
      idp_maximum(error + c, fabs(sum[c]) / rmax(1.0f, scale[c]));
  }
}

IDP_KERNEL void idp_floor_timestep_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                                          IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                                          IDP_GLOBAL const real *energy,
                                          IDP_GLOBAL const real *r0, IDP_GLOBAL const real *r1,
                                          IDP_GLOBAL const real *r2, IDP_GLOBAL const real *r3,
                                          IDP_GLOBAL const real *r4, IDP_GLOBAL real *limit,
                                          IDP_VALUE(real, floor), IDP_VALUE(real, upper),
                                          IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    real s[5], r[5];
    load_state(rho, mx, my, mz, energy, i, s);
    load_state(r0, r1, r2, r3, r4, i, r);
    limit[i] = floor_timestep(s, r, floor, upper);
  }
}

IDP_KERNEL void idp_low_update_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                                      IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                                      IDP_GLOBAL const real *energy, IDP_GLOBAL real *r0,
                                      IDP_GLOBAL real *r1, IDP_GLOBAL real *r2,
                                      IDP_GLOBAL real *r3, IDP_GLOBAL real *r4,
                                      IDP_VALUE(real, dt), IDP_VALUE(int, scalar_mode),
                                      IDP_VALUE(int, n) IDP_GRID) {
  IDP_GLOBAL real *out[5] = {r0, r1, r2, r3, r4};
  IDP_GLOBAL const real *in[5] = {rho, mx, my, mz, energy};
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS)
    for (int c = 0; c < 5; ++c)
      out[c][i] = (scalar_mode && c > 0) ? in[c][i] : in[c][i] - dt * out[c][i];
}

IDP_KERNEL void idp_bounds_init_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                                       IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                                       IDP_GLOBAL const real *energy, IDP_GLOBAL real *lower,
                                       IDP_GLOBAL real *upper, IDP_GLOBAL real *entropy,
                                       IDP_GLOBAL real *work, IDP_VALUE(real, gamma),
                                       IDP_VALUE(int, use_entropy),
                                       IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    lower[i] = upper[i] = rho[i];
    work[i] = 0.0f;
    if (use_entropy) {
      real s[5];
      load_state(rho, mx, my, mz, energy, i, s);
      entropy[i] = specific_entropy(s, gamma);
    }
  }
}

IDP_KERNEL void idp_bounds_edge_kernel(
    IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx, IDP_GLOBAL const real *my,
    IDP_GLOBAL const real *mz, IDP_GLOBAL const int *left, IDP_GLOBAL const int *right,
    IDP_GLOBAL const real *coef, IDP_GLOBAL const real *visc,
    IDP_GLOBAL const real *second_weight_left, IDP_GLOBAL const real *second_weight_right,
    IDP_GLOBAL real *lower, IDP_GLOBAL real *upper, IDP_GLOBAL const real *entropy_stage,
    IDP_GLOBAL real *entropy, IDP_GLOBAL real *work, IDP_VALUE(int, relax),
    IDP_VALUE(int, use_entropy), IDP_VALUE(int, nedge) IDP_GRID) {
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nedge; e += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[e], b = right[e];
    const real ra = rho[a], rb = rho[b];
    const real fa = coef[3 * e] * mx[a] + coef[3 * e + 1] * my[a] + coef[3 * e + 2] * mz[a];
    const real fb = coef[3 * e] * mx[b] + coef[3 * e + 1] * my[b] + coef[3 * e + 2] * mz[b];
    real bar = 0.5f * (ra + rb);
    if (visc[e] > IDP_TINY)
      bar -= (fb - fa) / (2.0f * visc[e]);
    idp_minimum(lower + a, rmin(rb, bar));
    idp_minimum(lower + b, rmin(ra, bar));
    idp_maximum(upper + a, rmax(rb, bar));
    idp_maximum(upper + b, rmax(ra, bar));
    if (use_entropy) {
      const real ea = entropy_stage[a], eb = entropy_stage[b];
      idp_minimum(entropy + a, eb);
      idp_minimum(entropy + b, ea);
    }
    if (relax) {
      const real diff = ra - rb;
      idp_add(work + a, second_weight_left[e] * diff);
      idp_add(work + b, -second_weight_right[e] * diff);
    }
  }
}

IDP_KERNEL void idp_relax_edge_kernel(IDP_GLOBAL const int *left, IDP_GLOBAL const int *right,
                                      IDP_GLOBAL const int *direction,
                                      IDP_GLOBAL const real *degree0,
                                      IDP_GLOBAL const real *degree1,
                                      IDP_GLOBAL const real *degree2,
                                      IDP_GLOBAL const real *first, IDP_GLOBAL real *second,
                                      IDP_VALUE(int, nedge) IDP_GRID) {
  IDP_GLOBAL const real *degree[3] = {degree0, degree1, degree2};
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nedge; e += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[e], b = right[e], d = direction[e] - 1;
    const real avg = 0.5f * (first[a] + first[b]);
    idp_add(second + a, 2.0f / degree[d][a] * avg);
    idp_add(second + b, 2.0f / degree[d][b] * avg);
  }
}

IDP_KERNEL void idp_relax_finalize_kernel(IDP_GLOBAL real *lower, IDP_GLOBAL real *upper,
                                          IDP_GLOBAL const real *second,
                                          IDP_VALUE(real, factor), IDP_VALUE(real, nodal_mass),
                                          IDP_VALUE(real, volume), IDP_VALUE(int, dimensions),
                                          IDP_VALUE(int, n) IDP_GRID) {
  const real rh = pow(nodal_mass / volume, 1.5f / (real)dimensions);
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    const real d = second[i] / (2.0f * (real)(2 * dimensions + 1));
    const real relax = factor * fabs(d);
    lower[i] = rmax((1.0f - rh) * lower[i], lower[i] - relax);
    upper[i] += relax;
  }
}

IDP_KERNEL void
idp_entropy_relax_edge_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                              IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                              IDP_GLOBAL const real *energy, IDP_GLOBAL const int *left,
                              IDP_GLOBAL const int *right, IDP_GLOBAL real *midpoint_max,
                              IDP_VALUE(real, gamma), IDP_VALUE(int, nedge) IDP_GRID) {
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nedge; e += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[e], b = right[e];
    real sa[5], sb[5], midpoint[5];
    load_state(rho, mx, my, mz, energy, a, sa);
    load_state(rho, mx, my, mz, energy, b, sb);
    for (int c = 0; c < 5; ++c)
      midpoint[c] = 0.5f * (sa[c] + sb[c]);
    const real entropy = specific_entropy(midpoint, gamma);
    idp_maximum(midpoint_max + a, entropy);
    idp_maximum(midpoint_max + b, entropy);
  }
}

IDP_KERNEL void
idp_entropy_relax_finalize_kernel(IDP_GLOBAL real *lower, IDP_GLOBAL const real *midpoint_max,
                                  IDP_VALUE(real, factor), IDP_VALUE(real, nodal_mass),
                                  IDP_VALUE(real, volume), IDP_VALUE(int, dimensions),
                                  IDP_VALUE(int, n) IDP_GRID) {
  const real rh = pow(nodal_mass / volume, 1.5f / (real)dimensions);
  const real cutoff = log(1.0f + rh / factor);
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    const real difference = rmax(0.0f, midpoint_max[i] - lower[i]);
    const real variation = exp(rmin(difference, cutoff)) - 1.0f;
    const real relaxation = rmin(rh, factor * rmax(0.0f, variation));
    lower[i] += log(1.0f - relaxation);
  }
}

IDP_KERNEL void idp_blend_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx,
                                 IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                                 IDP_GLOBAL const real *energy, IDP_GLOBAL const int *left,
                                 IDP_GLOBAL const int *right, IDP_GLOBAL const real *visc,
                                 IDP_GLOBAL const real *ev, IDP_GLOBAL real *correction,
                                 IDP_VALUE(real, dt), IDP_VALUE(int, has_ev),
                                 IDP_VALUE(int, low_only), IDP_VALUE(int, scalar_mode),
                                 IDP_VALUE(int, nedge) IDP_GRID) {
  IDP_GLOBAL const real *u[5] = {rho, mx, my, mz, energy};
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nedge; e += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[e], b = right[e];
    real high = low_only ? 0.0f : 1.0f;
    if (!low_only && has_ev)
      high = 1.0f - rmin(1.0f, rmax(0.0f, rmax(ev[a], ev[b])));
    for (int c = 0; c < 5; ++c)
      correction[5 * e + c] =
          (scalar_mode && c > 0)
              ? 0.0f
              : dt * high * (correction[5 * e + c] + visc[e] * (u[c][a] - u[c][b]));
  }
}

IDP_KERNEL void idp_limiter_kernel(
    IDP_GLOBAL const real *q0, IDP_GLOBAL const real *q1, IDP_GLOBAL const real *q2,
    IDP_GLOBAL const real *q3, IDP_GLOBAL const real *q4, IDP_GLOBAL const real *lower,
    IDP_GLOBAL const real *upper, IDP_GLOBAL const real *entropy, IDP_GLOBAL const real *mass,
    IDP_GLOBAL const real *degree0, IDP_GLOBAL const real *degree1,
    IDP_GLOBAL const real *degree2, IDP_GLOBAL const int *left, IDP_GLOBAL const int *right,
    IDP_GLOBAL const int *direction, IDP_GLOBAL real *correction, IDP_GLOBAL real *edge_limit,
    IDP_GLOBAL real *limited, IDP_GLOBAL real *density_flag, IDP_GLOBAL real *energy_flag,
    IDP_GLOBAL real *entropy_flag, IDP_VALUE(real, gamma), IDP_VALUE(real, floor),
    IDP_VALUE(int, enforce_energy), IDP_VALUE(int, enforce_entropy), IDP_VALUE(int, check_base),
    IDP_VALUE(int, diagnostics_level), IDP_VALUE(int, dimensions),
    IDP_VALUE(int, nedge) IDP_GRID) {
  IDP_GLOBAL const real *q[5] = {q0, q1, q2, q3, q4};
  IDP_GLOBAL const real *degree[3] = {degree0, degree1, degree2};
  const real eps = IDP_EPS;
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nedge; e += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[e], b = right[e], d = direction[e] - 1;
    real lb[5], rb[5], lc[5], rc[5];
    real scale = 1.0f, size = 0.0f;
    for (int c = 0; c < 5; ++c) {
      lb[c] = q[c][a];
      rb[c] = q[c][b];
      lc[c] = correction[5 * e + c] * (real)dimensions * degree[d][a] / mass[a];
      rc[c] = -correction[5 * e + c] * (real)dimensions * degree[d][b] / mass[b];
      scale = rmax(scale, rmax(fabs(lb[c]), fabs(rb[c])));
      size = rmax(size, rmax(fabs(lc[c]), fabs(rc[c])));
    }
    bool dl = false, el = false, sl = false;
    real limit = 1.0f;
    if (size <= 512.0f * eps * scale) {
      for (int c = 0; c < 5; ++c)
        correction[5 * e + c] = 0.0f;
    } else {
      bool dl1, el1, sl1, dl2, el2, sl2;
      const real sent = enforce_entropy ? entropy[a] : -IDP_HUGE;
      const real tent = enforce_entropy ? entropy[b] : -IDP_HUGE;
      real l1 = limit_endpoint(lb, lc, lower[a], upper[a], sent, gamma, floor, enforce_energy,
                               enforce_entropy, check_base, &dl1, &el1, &sl1);
      real l2 = limit_endpoint(rb, rc, lower[b], upper[b], tent, gamma, floor, enforce_energy,
                               enforce_entropy, check_base, &dl2, &el2, &sl2);
      limit = rmin(l1, l2);
      dl = dl1 || dl2;
      el = el1 || el2;
      sl = sl1 || sl2;
      if (enforce_entropy && limit > 0.0f && limit < 1.0f) {
        real lt[5], rt[5];
        for (int c = 0; c < 5; ++c) {
          lt[c] = lb[c] + limit * lc[c];
          rt[c] = rb[c] + limit * rc[c];
        }
        if (!entropy_is_admissible(lt, gamma, sent, floor) ||
            !entropy_is_admissible(rt, gamma, tent, floor)) {
          real lo = 0.0f, hi = limit;
          for (int it = 0; it < 64; ++it) {
            const real mid = lo + 0.5f * (hi - lo);
            for (int c = 0; c < 5; ++c) {
              lt[c] = lb[c] + mid * lc[c];
              rt[c] = rb[c] + mid * rc[c];
            }
            if (entropy_is_admissible(lt, gamma, sent, floor) &&
                entropy_is_admissible(rt, gamma, tent, floor))
              lo = mid;
            else
              hi = mid;
          }
          limit = lo;
          sl = true;
        }
      }
      for (int c = 0; c < 5; ++c)
        correction[5 * e + c] *= limit;
    }
    if (diagnostics_level > 0)
      edge_limit[e] = limit;
    if (diagnostics_level > 1) {
      limited[e] = limit < 1.0f - 32.0f * eps ? 1.0f : 0.0f;
      density_flag[e] = dl ? 1.0f : 0.0f;
      energy_flag[e] = el ? 1.0f : 0.0f;
      entropy_flag[e] = sl ? 1.0f : 0.0f;
    }
  }
}

IDP_KERNEL void idp_incidence_kernel(IDP_GLOBAL const int *left, IDP_GLOBAL const int *right,
                                     IDP_GLOBAL const real *correction, IDP_GLOBAL real *r0,
                                     IDP_GLOBAL real *r1, IDP_GLOBAL real *r2,
                                     IDP_GLOBAL real *r3, IDP_GLOBAL real *r4,
                                     IDP_VALUE(int, nedge) IDP_GRID) {
  IDP_GLOBAL real *r[5] = {r0, r1, r2, r3, r4};
  for (int e = IDP_GROUP * IDP_WIDTH + IDP_LANE; e < nedge; e += IDP_WIDTH * IDP_GROUPS) {
    const int a = left[e], b = right[e];
    for (int c = 0; c < 5; ++c) {
      idp_add(r[c] + a, correction[5 * e + c]);
      idp_add(r[c] + b, -correction[5 * e + c]);
    }
  }
}

IDP_KERNEL void idp_correction_update_kernel(
    IDP_GLOBAL real *q0, IDP_GLOBAL real *q1, IDP_GLOBAL real *q2, IDP_GLOBAL real *q3,
    IDP_GLOBAL real *q4, IDP_GLOBAL const real *r0, IDP_GLOBAL const real *r1,
    IDP_GLOBAL const real *r2, IDP_GLOBAL const real *r3, IDP_GLOBAL const real *r4,
    IDP_GLOBAL const real *mass, IDP_VALUE(int, n) IDP_GRID) {
  IDP_GLOBAL real *q[5] = {q0, q1, q2, q3, q4};
  IDP_GLOBAL const real *r[5] = {r0, r1, r2, r3, r4};
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS)
    for (int c = 0; c < 5; ++c)
      q[c][i] += r[c][i] / mass[i];
}

IDP_KERNEL void idp_validation_kernel(IDP_GLOBAL const real *q0, IDP_GLOBAL const real *q1,
                                      IDP_GLOBAL const real *q2, IDP_GLOBAL const real *q3,
                                      IDP_GLOBAL const real *q4, IDP_GLOBAL const real *lower,
                                      IDP_GLOBAL const real *upper,
                                      IDP_GLOBAL const real *entropy_lower,
                                      IDP_GLOBAL real *lower_v, IDP_GLOBAL real *upper_v,
                                      IDP_GLOBAL real *entropy_v,
                                      IDP_GLOBAL real *entropy_value, IDP_VALUE(real, gamma),
                                      IDP_VALUE(int, use_entropy), IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    real s[5];
    load_state(q0, q1, q2, q3, q4, i, s);
    lower_v[i] = rmax(0.0f, lower[i] - q0[i]);
    upper_v[i] = rmax(0.0f, q0[i] - upper[i]);
    const real ent = specific_entropy(s, gamma);
    entropy_value[i] = ent;
    entropy_v[i] =
        use_entropy
            ? rmax(0.0f, entropy_lower[i] - ent - entropy_tolerance(s, entropy_lower[i]))
            : 0.0f;
  }
}

enum diagnostic_summary_kind {
  validation_summary = 1,
  observation_summary = 2,
  limiter_status_summary = 3
};

IDP_KERNEL void idp_diagnostic_summary_init_kernel(IDP_GLOBAL real *summary,
                                                   IDP_VALUE(int, size),
                                                   IDP_VALUE(int, kind) IDP_GRID) {
  if (IDP_GROUP != 0 || IDP_LANE != 0)
    return;
  for (int i = 0; i < size; ++i)
    summary[i] = 0.0f;
  if (kind == validation_summary) {
    for (int i = 3; i <= 6; ++i)
      summary[i] = IDP_HUGE;
  } else if (kind == observation_summary) {
    for (int i = 0; i <= 2; ++i)
      summary[i] = IDP_HUGE;
  } else if (kind == limiter_status_summary) {
    summary[0] = IDP_HUGE;
    summary[1] = -IDP_HUGE;
  }
}

enum full_diagnostic_partial_size { full_node_partial_size = 8, full_edge_partial_size = 10 };

IDP_KERNEL void idp_full_diagnostics_node_partial_kernel(
    IDP_GLOBAL const real *rho, IDP_GLOBAL const real *mx, IDP_GLOBAL const real *my,
    IDP_GLOBAL const real *mz, IDP_GLOBAL const real *energy, IDP_VALUE(real, gamma),
    IDP_GLOBAL const real *entropy_fraction, IDP_VALUE(int, has_entropy),
    IDP_GLOBAL const real *r0, IDP_GLOBAL const real *r1, IDP_GLOBAL const real *r2,
    IDP_GLOBAL const real *r3, IDP_GLOBAL const real *r4, IDP_GLOBAL real *partial,
    IDP_VALUE(int, n) IDP_GRID) {
  IDP_SHARED real shared[256 * full_node_partial_size];

  IDP_GLOBAL const real *residual[5] = {r0, r1, r2, r3, r4};
  const int tid = IDP_LANE;
  const int i = IDP_GROUP * IDP_WIDTH + tid;
  real value[full_node_partial_size] = {0.0f, 0.0f, 0.0f, 0.0f,
                                        0.0f, 0.0f, 0.0f, 0.0f};
  if (i < n) {
    real state[5];
    load_state(rho, mx, my, mz, energy, i, state);
    const real pressure = (gamma - 1.0f) * internal_energy(state);
    const real velocity = sqrt(mx[i] * mx[i] + my[i] * my[i] + mz[i] * mz[i]) / rho[i];
    value[0] = velocity + sqrt(gamma * pressure / rho[i]);
    if (has_entropy) {
      value[1] = entropy_fraction[i];
      value[2] = entropy_fraction[i];
    }
    for (int c = 0; c < 5; ++c)
      value[3 + c] = residual[c][i];
  }
  for (int metric = 0; metric < full_node_partial_size; ++metric)
    shared[metric * IDP_WIDTH + tid] = value[metric];
  IDP_BARRIER();
  for (int stride = IDP_WIDTH / 2; stride > 0; stride /= 2) {
    if (tid < stride) {
      shared[tid] = rmax(shared[tid], shared[tid + stride]);
      shared[IDP_WIDTH + tid] = rmax(shared[IDP_WIDTH + tid], shared[IDP_WIDTH + tid + stride]);
      for (int metric = 2; metric < full_node_partial_size; ++metric)
        shared[metric * IDP_WIDTH + tid] += shared[metric * IDP_WIDTH + tid + stride];
    }
    IDP_BARRIER();
  }
  if (tid < full_node_partial_size)
    partial[IDP_GROUP * full_node_partial_size + tid] = shared[tid * IDP_WIDTH];
}

IDP_KERNEL void idp_full_diagnostics_edge_partial_kernel(
    IDP_GLOBAL const real *edge_limit, IDP_GLOBAL const real *limited,
    IDP_GLOBAL const real *density_flag, IDP_GLOBAL const real *energy_flag,
    IDP_GLOBAL const real *entropy_flag, IDP_GLOBAL const real *correction,
    IDP_GLOBAL real *partial, IDP_VALUE(int, nedge), IDP_VALUE(int, node_blocks) IDP_GRID) {
  IDP_SHARED real shared[256 * full_edge_partial_size];

  const int tid = IDP_LANE;
  const int i = IDP_GROUP * IDP_WIDTH + tid;
  const real huge = IDP_HUGE;
  const real eps = IDP_EPS;
  real value[full_edge_partial_size] = {huge,        -huge,       0.0f, 0.0f,
                                        0.0f, 0.0f, 0.0f, 0.0f,
                                        0.0f, 0.0f};
  if (i < nedge) {
    const real limit = edge_limit[i];
    value[0] = limit;
    value[1] = limit;
    value[2] = limit;
    value[3] = limited[i];
    value[4] = density_flag[i];
    value[5] = energy_flag[i];
    value[6] = entropy_flag[i];
    value[7] = (!isfinite(limit) || limit < -32.0f * eps ||
                limit > 1.0f + 32.0f * eps)
                   ? 1.0f
                   : 0.0f;
    for (int c = 0; c < 5; ++c) {
      const real flux = correction[5 * i + c];
      value[8] = rmax(value[8], fabs(flux));
      value[9] += flux * flux;
    }
  }
  for (int metric = 0; metric < full_edge_partial_size; ++metric)
    shared[metric * IDP_WIDTH + tid] = value[metric];
  IDP_BARRIER();
  for (int stride = IDP_WIDTH / 2; stride > 0; stride /= 2) {
    if (tid < stride) {
      shared[tid] = rmin(shared[tid], shared[tid + stride]);
      shared[IDP_WIDTH + tid] = rmax(shared[IDP_WIDTH + tid], shared[IDP_WIDTH + tid + stride]);
      for (int metric = 2; metric < 8; ++metric)
        shared[metric * IDP_WIDTH + tid] += shared[metric * IDP_WIDTH + tid + stride];
      shared[8 * IDP_WIDTH + tid] =
          rmax(shared[8 * IDP_WIDTH + tid], shared[8 * IDP_WIDTH + tid + stride]);
      shared[9 * IDP_WIDTH + tid] += shared[9 * IDP_WIDTH + tid + stride];
    }
    IDP_BARRIER();
  }
  if (tid < full_edge_partial_size)
    partial[full_node_partial_size * node_blocks + IDP_GROUP * full_edge_partial_size + tid] =
        shared[tid * IDP_WIDTH];
}

IDP_KERNEL void idp_full_diagnostics_finalize_kernel(IDP_GLOBAL const real *node_partial,
                                                     IDP_GLOBAL const real *directional_error,
                                                     IDP_GLOBAL real *summary,
                                                     IDP_VALUE(int, node_blocks),
                                                     IDP_VALUE(int, edge_blocks) IDP_GRID) {
  const int metric = IDP_LANE;
  if (IDP_GROUP != 0 || metric >= 23)
    return;
  real result = 0.0f;
  if (metric == 0 || metric == 1) {
    for (int block = 0; block < node_blocks; ++block)
      result = rmax(result, node_partial[block * full_node_partial_size + metric]);
  } else if (metric == 2) {
    for (int block = 0; block < node_blocks; ++block)
      result += node_partial[block * full_node_partial_size + 2];
  } else if (metric >= 10 && metric < 15) {
    const int source = metric - 7;
    for (int block = 0; block < node_blocks; ++block)
      result += node_partial[block * full_node_partial_size + source];
  } else if (metric >= 15 && metric < 20) {
    result = directional_error[metric - 15];
  } else {
    int source = -1;
    if (metric >= 3 && metric <= 9)
      source = metric - 3;
    else if (metric == 20)
      source = 8;
    else if (metric == 21)
      source = 9;
    else if (metric == 22)
      source = 7;
    if (metric == 3)
      result = IDP_HUGE;
    else if (metric == 4)
      result = -IDP_HUGE;
    IDP_GLOBAL const real *edge_partial = node_partial + full_node_partial_size * node_blocks;
    for (int block = 0; block < edge_blocks; ++block) {
      const real value = edge_partial[block * full_edge_partial_size + source];
      if (metric == 3)
        result = rmin(result, value);
      else if (metric == 4 || metric == 20)
        result = rmax(result, value);
      else
        result += value;
    }
  }
  summary[metric] = result;
}

IDP_KERNEL void idp_validation_summary_kernel(
    IDP_GLOBAL const real *q0, IDP_GLOBAL const real *q1, IDP_GLOBAL const real *q2,
    IDP_GLOBAL const real *q3, IDP_GLOBAL const real *q4, IDP_GLOBAL const real *lower,
    IDP_GLOBAL const real *upper, IDP_GLOBAL const real *entropy_lower,
    IDP_GLOBAL real *summary, IDP_VALUE(real, gamma), IDP_VALUE(int, use_entropy),
    IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    real state[5];
    load_state(q0, q1, q2, q3, q4, i, state);
    const real internal = internal_energy(state);
    const real pressure = (gamma - 1.0f) * internal;
    const real entropy = specific_entropy(state, gamma);
    const real lower_violation = rmax(0.0f, lower[i] - q0[i]);
    const real upper_violation = rmax(0.0f, q0[i] - upper[i]);
    const real entropy_violation =
        use_entropy ? rmax(0.0f, entropy_lower[i] - entropy -
                                            entropy_tolerance(state, entropy_lower[i]))
                    : 0.0f;
    idp_maximum(summary, lower_violation);
    idp_maximum(summary + 1, upper_violation);
    idp_maximum(summary + 2, entropy_violation);
    idp_minimum(summary + 3, q0[i]);
    idp_minimum(summary + 4, internal);
    idp_minimum(summary + 5, pressure);
    idp_minimum(summary + 6, entropy);
    idp_maximum(summary + 7, rmax(q0[i], upper[i]));
    bool invalid = q0[i] <= 0.0f || !isfinite(internal) || internal <= 0.0f ||
                   !isfinite(pressure) || pressure <= 0.0f || !isfinite(entropy);
    for (int c = 0; c < 5; ++c)
      invalid = invalid || !isfinite(state[c]);
    if (invalid)
      idp_add(summary + 8, 1.0f);
  }
}

IDP_KERNEL void
idp_observation_summary_kernel(IDP_GLOBAL const real *rho, IDP_GLOBAL const real *internal,
                               IDP_GLOBAL const real *pressure, IDP_GLOBAL const real *u,
                               IDP_GLOBAL const real *v, IDP_GLOBAL const real *w,
                               IDP_GLOBAL const real *sound, IDP_GLOBAL real *summary,
                               IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    idp_minimum(summary, rho[i]);
    idp_minimum(summary + 1, internal[i]);
    idp_minimum(summary + 2, pressure[i]);
    const real wave = sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]) + sound[i];
    idp_maximum(summary + 3, wave);
  }
}

IDP_KERNEL void idp_limiter_status_kernel(IDP_GLOBAL const real *edge_limit,
                                          IDP_GLOBAL real *summary,
                                          IDP_VALUE(int, nedge) IDP_GRID) {
  const real eps = IDP_EPS;
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < nedge; i += IDP_WIDTH * IDP_GROUPS) {
    const real limit = edge_limit[i];
    idp_minimum(summary, limit);
    idp_maximum(summary + 1, limit);
    if (!isfinite(limit) || limit < -32.0f * eps ||
        limit > 1.0f + 32.0f * eps)
      idp_add(summary + 2, 1.0f);
  }
}

IDP_KERNEL void idp_update_uvw_kernel_idp(IDP_GLOBAL real *u, IDP_GLOBAL real *v,
                                          IDP_GLOBAL real *w, IDP_GLOBAL const real *mx,
                                          IDP_GLOBAL const real *my, IDP_GLOBAL const real *mz,
                                          IDP_GLOBAL const real *rho,
                                          IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    u[i] = mx[i] / rho[i];
    v[i] = my[i] / rho[i];
    w[i] = mz[i] / rho[i];
  }
}

IDP_KERNEL void idp_update_momentum_pressure_kernel(
    IDP_GLOBAL real *mx, IDP_GLOBAL real *my, IDP_GLOBAL real *mz, IDP_GLOBAL real *p,
    IDP_GLOBAL real *kin, IDP_GLOBAL const real *u, IDP_GLOBAL const real *v,
    IDP_GLOBAL const real *w, IDP_GLOBAL const real *energy, IDP_GLOBAL const real *rho,
    IDP_VALUE(real, gamma), IDP_VALUE(int, n) IDP_GRID) {
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    mx[i] = rho[i] * u[i];
    my[i] = rho[i] * v[i];
    mz[i] = rho[i] * w[i];
    kin[i] = 0.5f * rho[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    p[i] = (gamma - 1.0f) * (energy[i] - kin[i]);
  }
}

IDP_KERNEL void idp_update_energy_kernel(IDP_GLOBAL real *energy, IDP_GLOBAL real *p,
                                         IDP_GLOBAL const real *kin, IDP_VALUE(real, gamma),
                                         IDP_VALUE(real, input_floor),
                                         IDP_VALUE(int, n) IDP_GRID) {
  const real inv = 1.0f / (gamma - 1.0f);
  const real floor = rmax(input_floor, 0.0f);
  for (int i = IDP_GROUP * IDP_WIDTH + IDP_LANE; i < n; i += IDP_WIDTH * IDP_GROUPS) {
    if (p[i] * inv <= floor) {
      p[i] = rmax(p[i], (real)1e-12);
      if (p[i] * inv <= floor + IDP_EPS * fabs(kin[i])) {
        const real guard =
            32.0f * IDP_EPS * rmax(1.0f, rmax(fabs(kin[i]), floor));
        p[i] = (gamma - 1.0f) * (floor + guard);
      }
    }
    energy[i] = p[i] * inv + kin[i];
  }
}
