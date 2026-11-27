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

/*
 * Trajectory-based collision of particles with the high-order wall facets,
 * shared by the CUDA and HIP backends. This is a line-by-line port of
 * lpt_wall_geometry.f90 (see there for the formulation); every function
 * has the Fortran routine of the same name as its reference.
 *
 * No dynamic memory is used; all work arrays have compile-time sizes, so
 * lx must not exceed LPT_WALL_MAX_LX.
 */

#ifndef __LPT_WALL_INTERSECT_H__
#define __LPT_WALL_INTERSECT_H__

#include <math.h>

#define LPT_WALL_MAX_LX 16
#define LPT_WALL_MAX_BOUNCES 8
#define LPT_NEWTON_MAX_ITER 25
#define LPT_NEWTON_MAX_HALVINGS 10
#define LPT_PATCH_MAX_DEPTH 4
#define LPT_PATCH_STACK (3 * LPT_PATCH_MAX_DEPTH + 2)

/* Newton outcomes */
#define LPT_NEWTON_CONVERGED 0
#define LPT_NEWTON_SINGULAR 1
#define LPT_NEWTON_NOT_CONVERGED 2
#define LPT_NEWTON_OUTSIDE 3

/* Layout of the per-facet real record, see lpt_wall_geometry_t%to_device */
#define LPT_FR_ORIGIN 0
#define LPT_FR_ORIENT 3
#define LPT_FR_CENTER 4   /* facet-local */
#define LPT_FR_RADIUS 7
#define LPT_FR_AXIS 8
#define LPT_FR_ANGLE 11
#define LPT_FR_SLAB 12
#define LPT_FR_SCALE 14
#define LPT_FR_PLANAR 15
#define LPT_FR_ELEN 16
#define LPT_FR_EC 20      /* 3 x 4, facet-local */
#define LPT_FR_ER 32
#define LPT_FR_EB 36
#define LPT_FR_CX 40      /* 3 x 4, facet-local */
#define LPT_FR_CN 52      /* 3 x 4 */
#define LPT_NFR 64

/* Diagnostic counters, same order as lpt_wall_stats_t on the host */
#define LPT_ST_UNRESOLVED 0
#define LPT_ST_REACH 1
#define LPT_NSTATS 2

template<typename T> __device__ inline T lpt_eps();
template<> __device__ inline float lpt_eps<float>() { return 1.1920929e-7f; }
template<> __device__ inline double lpt_eps<double>() { return 2.220446049250313e-16; }
template<typename T> __device__ inline T lpt_huge();
template<> __device__ inline float lpt_huge<float>() { return 1.0e30f; }
template<> __device__ inline double lpt_huge<double>() { return 1.0e300; }

template<typename T>
struct lpt_v3 {
  T x, y, z;
};

template<typename T>
__device__ inline lpt_v3<T> mk3(T x, T y, T z) { lpt_v3<T> r; r.x = x; r.y = y; r.z = z; return r; }
template<typename T>
__device__ inline lpt_v3<T> operator+(lpt_v3<T> a, lpt_v3<T> b) { return mk3(a.x + b.x, a.y + b.y, a.z + b.z); }
template<typename T>
__device__ inline lpt_v3<T> operator-(lpt_v3<T> a, lpt_v3<T> b) { return mk3(a.x - b.x, a.y - b.y, a.z - b.z); }
template<typename T>
__device__ inline lpt_v3<T> operator*(T s, lpt_v3<T> a) { return mk3(s * a.x, s * a.y, s * a.z); }
template<typename T>
__device__ inline T dot3(lpt_v3<T> a, lpt_v3<T> b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
template<typename T>
__device__ inline T nrm3(lpt_v3<T> a) { return sqrt(dot3(a, a)); }
template<typename T>
__device__ inline lpt_v3<T> cross3(lpt_v3<T> a, lpt_v3<T> b) {
  return mk3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x);
}
template<typename T>
__device__ inline T lpt_max(T a, T b) { return a > b ? a : b; }
template<typename T>
__device__ inline T lpt_min(T a, T b) { return a < b ? a : b; }
template<typename T>
__device__ inline T lpt_clamp(T v, T lo, T hi) { return lpt_max(lo, lpt_min(hi, v)); }

/* Read-only geometry, packed by lpt_wall_geometry_t%to_device */
template<typename T>
struct lpt_wall_geom {
  const T *xf;          /* (3, lx, lx, nf), facet-local */
  const T *fr;          /* (LPT_NFR, nf) */
  const int *fi;        /* (2, nf): local element, element facet */
  const T *zg;          /* (lx) */
  const int *el_ptr;    /* (nelv + 1), 1-based CSR */
  const int *el_facets; /* 1-based facet indices */
  const T *el_radius;   /* (nelv) */
  int lx;
  int nelv;
};

/* Per-thread counters, added to the global ones once per particle */
struct lpt_wall_cnt {
  unsigned int c[LPT_NSTATS];
};

template<typename T>
__device__ inline const T *fr_of(const lpt_wall_geom<T> &g, int f) {
  return g.fr + (size_t) LPT_NFR * (f - 1);
}

template<typename T>
__device__ inline lpt_v3<T> fr3(const lpt_wall_geom<T> &g, int f, int off) {
  const T *r = fr_of(g, f);
  return mk3(r[off], r[off + 1], r[off + 2]);
}

/* Lagrange weights and derivatives up to order m at xi (Fornberg), as
   fd_weights_full; computed in double. c[j][k], j node, k derivative. */
template<typename T>
__device__ inline void lpt_fd_weights(T xi_, const T *zg, int lx, int m,
                            double c[LPT_WALL_MAX_LX][3]) {
  const double xi = (double) xi_;
  const int n = lx - 1;
  double c1 = 1.0, c4 = (double) zg[0] - xi;
  for (int j = 0; j <= n; j++)
    for (int k = 0; k <= m; k++)
      c[j][k] = 0.0;
  c[0][0] = 1.0;
  for (int i = 1; i <= n; i++) {
    const int mn = i < m ? i : m;
    double c2 = 1.0;
    const double c5 = c4;
    c4 = (double) zg[i] - xi;
    for (int j = 0; j <= i - 1; j++) {
      const double c3 = (double) zg[i] - (double) zg[j];
      c2 = c2 * c3;
      for (int k = mn; k >= 1; k--)
        c[i][k] = c1 * (k * c[i - 1][k - 1] - c5 * c[i - 1][k]) / c2;
      c[i][0] = -c1 * c5 * c[i - 1][0] / c2;
      for (int k = mn; k >= 1; k--)
        c[j][k] = (c4 * c[j][k] - k * c[j][k - 1]) / c3;
      c[j][0] = c4 * c[j][0] / c3;
    }
    c1 = c2;
  }
}

template<typename T>
__device__ inline lpt_v3<T> xf_node(const lpt_wall_geom<T> &g, int f, int a, int b) {
  /* a, b 1-based */
  const size_t i = 3 * ((size_t) (a - 1) + (size_t) g.lx *
                        ((size_t) (b - 1) + (size_t) g.lx * (f - 1)));
  return mk3(g.xf[i], g.xf[i + 1], g.xf[i + 2]);
}

/* eval_local / eval_nodes */
template<typename T>
__device__ inline void lpt_eval2(const lpt_wall_geom<T> &g, int f, T s, T t,
                       lpt_v3<T> &x, lpt_v3<T> &xs, lpt_v3<T> &xt,
                       lpt_v3<T> &xss, lpt_v3<T> &xst, lpt_v3<T> &xtt) {
  double cs[LPT_WALL_MAX_LX][3], ct[LPT_WALL_MAX_LX][3];
  lpt_fd_weights<T>(s, g.zg, g.lx, 2, cs);
  lpt_fd_weights<T>(t, g.zg, g.lx, 2, ct);
  const lpt_v3<T> z = mk3<T>(0, 0, 0);
  x = z; xs = z; xt = z; xss = z; xst = z; xtt = z;
  for (int b = 1; b <= g.lx; b++) {
    lpt_v3<T> a0 = z, a1 = z, a2 = z;
    for (int a = 1; a <= g.lx; a++) {
      const lpt_v3<T> p = xf_node(g, f, a, b);
      a0 = a0 + (T) cs[a - 1][0] * p;
      a1 = a1 + (T) cs[a - 1][1] * p;
      a2 = a2 + (T) cs[a - 1][2] * p;
    }
    x = x + (T) ct[b - 1][0] * a0;
    xs = xs + (T) ct[b - 1][0] * a1;
    xt = xt + (T) ct[b - 1][1] * a0;
    xss = xss + (T) ct[b - 1][0] * a2;
    xst = xst + (T) ct[b - 1][1] * a1;
    xtt = xtt + (T) ct[b - 1][2] * a0;
  }
}

/* eval_first */
template<typename T>
__device__ inline void lpt_eval1(const lpt_wall_geom<T> &g, int f, T s, T t,
                       lpt_v3<T> &x, lpt_v3<T> &xs, lpt_v3<T> &xt) {
  double cs[LPT_WALL_MAX_LX][3], ct[LPT_WALL_MAX_LX][3];
  lpt_fd_weights<T>(s, g.zg, g.lx, 1, cs);
  lpt_fd_weights<T>(t, g.zg, g.lx, 1, ct);
  const lpt_v3<T> z = mk3<T>(0, 0, 0);
  x = z; xs = z; xt = z;
  for (int b = 1; b <= g.lx; b++) {
    lpt_v3<T> a0 = z, a1 = z;
    for (int a = 1; a <= g.lx; a++) {
      const lpt_v3<T> p = xf_node(g, f, a, b);
      a0 = a0 + (T) cs[a - 1][0] * p;
      a1 = a1 + (T) cs[a - 1][1] * p;
    }
    x = x + (T) ct[b - 1][0] * a0;
    xs = xs + (T) ct[b - 1][0] * a1;
    xt = xt + (T) ct[b - 1][1] * a0;
  }
}

/* surface_frame */
template<typename T>
__device__ inline bool lpt_frame(const lpt_wall_geom<T> &g, int f, T s, T t,
                       lpt_v3<T> &x, lpt_v3<T> &xs, lpt_v3<T> &xt,
                       lpt_v3<T> &n, bool deriv, lpt_v3<T> &ns,
                       lpt_v3<T> &nt) {
  lpt_v3<T> xss, xst, xtt;
  if (deriv)
    lpt_eval2(g, f, s, t, x, xs, xt, xss, xst, xtt);
  else
    lpt_eval1(g, f, s, t, x, xs, xt);
  const lpt_v3<T> nv = cross3(xs, xt);
  const T nlen = nrm3(nv);
  const T orient = fr_of(g, f)[LPT_FR_ORIENT];
  if (!(nlen > (T) 0)) {
    n = mk3<T>(0, 0, 0); ns = n; nt = n;
    return false;
  }
  const lpt_v3<T> nu = ((T) 1 / nlen) * nv;
  n = orient * nu;
  if (deriv) {
    const lpt_v3<T> ds = cross3(xss, xt) + cross3(xs, xst);
    const lpt_v3<T> dt = cross3(xst, xt) + cross3(xs, xtt);
    ns = (orient / nlen) * (ds - dot3(nu, ds) * nu);
    nt = (orient / nlen) * (dt - dot3(nu, dt) * nu);
  } else {
    ns = mk3<T>(0, 0, 0); nt = ns;
  }
  return true;
}

template<typename T>
__device__ inline T lpt_angle(lpt_v3<T> u, lpt_v3<T> v) {
  return atan2(nrm3(cross3(u, v)), dot3(u, v));
}

/* patch_bounds, in two passes over the samples (no sample storage) */
template<typename T>
__device__ inline void lpt_patch_bounds(const lpt_wall_geom<T> &g, int f,
                              const T box[4], lpt_v3<T> &c, T &r,
                              lpt_v3<T> &axis, T &theta, T slab[2]) {
  const int m = lpt_max(3, 2 * g.lx - 1);
  const T hs = (box[1] - box[0]) / (m - 1);
  const T ht = (box[3] - box[2]) / (m - 1);
  const T orient = fr_of(g, f)[LPT_FR_ORIENT];
  lpt_v3<T> sx = mk3<T>(0, 0, 0), sn = sx, x, xs, xt;
  T gs = 0, gt = 0;
  bool degenerate = false;
  for (int b = 1; b <= m; b++) {
    for (int a = 1; a <= m; a++) {
      lpt_eval1(g, f, box[0] + (a - 1) * hs, box[2] + (b - 1) * ht, x, xs,
                xt);
      gs = lpt_max(gs, nrm3(xs));
      gt = lpt_max(gt, nrm3(xt));
      const lpt_v3<T> nv = cross3(xs, xt);
      const T nl = nrm3(nv);
      if (nl > (T) 0) sn = sn + (orient / nl) * nv;
      else degenerate = true;
      sx = sx + x;
    }
  }
  c = ((T) 1 / (T) (m * m)) * sx;
  const T margin = (T) 0.5 * (hs * gs + ht * gt);
  const T nl = nrm3(sn);
  r = 0;
  if (degenerate || !(nl > (T) 0)) {
    for (int b = 1; b <= m; b++)
      for (int a = 1; a <= m; a++) {
        lpt_eval1(g, f, box[0] + (a - 1) * hs, box[2] + (b - 1) * ht, x,
                  xs, xt);
        r = lpt_max(r, nrm3(x - c));
      }
    r = (T) 1.01 * (r + margin);
    axis = mk3<T>(0, 0, 1);
    theta = lpt_huge<T>();
    slab[0] = -r;
    slab[1] = r;
    return;
  }
  axis = ((T) 1 / nl) * sn;
  theta = 0;
  T dmax = 0;
  slab[0] = lpt_huge<T>();
  slab[1] = -lpt_huge<T>();
  lpt_v3<T> prev_row[2 * LPT_WALL_MAX_LX + 1];
  for (int b = 1; b <= m; b++) {
    lpt_v3<T> left = mk3<T>(0, 0, 0);
    for (int a = 1; a <= m; a++) {
      lpt_eval1(g, f, box[0] + (a - 1) * hs, box[2] + (b - 1) * ht, x, xs,
                xt);
      r = lpt_max(r, nrm3(x - c));
      const T h = dot3(x - c, axis);
      slab[0] = lpt_min(slab[0], h);
      slab[1] = lpt_max(slab[1], h);
      const lpt_v3<T> nv = cross3(xs, xt);
      const lpt_v3<T> nn = (orient / nrm3(nv)) * nv;
      theta = lpt_max(theta, lpt_angle(nn, axis));
      if (a > 1) dmax = lpt_max(dmax, lpt_angle(nn, left));
      if (b > 1) dmax = lpt_max(dmax, lpt_angle(nn, prev_row[a - 1]));
      left = nn;
      prev_row[a - 1] = nn;
    }
  }
  r = (T) 1.01 * (r + margin);
  theta = (T) 1.01 * (theta + dmax);
  slab[0] -= (T) 1.01 * margin;
  slab[1] += (T) 1.01 * margin;
}

/* segment_sphere_interval */
template<typename T>
__device__ inline bool lpt_seg_sphere(lpt_v3<T> x0, lpt_v3<T> dx, lpt_v3<T> c, T r,
                            T &a_in, T &a_out) {
  a_in = 0; a_out = 0;
  const T qa = dot3(dx, dx);
  const T qb = dot3(dx, x0 - c);
  const T qc = dot3(x0 - c, x0 - c) - r * r;
  if (qa <= (T) 0) {
    a_out = 1;
    return qc <= (T) 0;
  }
  const T disc = qb * qb - qa * qc;
  if (disc < (T) 0) return false;
  const T sq = sqrt(disc);
  a_in = lpt_max((T) 0, (-qb - sq) / qa);
  a_out = lpt_min((T) 1, (-qb + sq) / qa);
  return a_in <= a_out;
}

template<typename T>
__device__ inline bool lpt_at_bound(T s, T t, T al, const T lo[2], const T hi[2]) {
  return s <= lo[0] || t <= lo[1] || s >= hi[0] || t >= hi[1] ||
    al <= (T) -0.5 || al >= (T) 1.5;
}

/* newton_intersect */
template<typename T>
__device__ inline int lpt_newton(const lpt_wall_geom<T> &g, int f, lpt_v3<T> x0,
                       lpt_v3<T> dx, T a, const T box[4], T s0, T t0,
                       T alpha0, bool exact, T &s, T &t, T &alpha,
                       lpt_v3<T> &n, int &iters) {
  const T eps = lpt_eps<T>();
  const T *rec = fr_of(g, f);
  const T ext = (T) 0.5 * lpt_max(box[1] - box[0], box[3] - box[2]);
  const T lo[2] = {lpt_max((T) -1.25, box[0] - ext),
                   lpt_max((T) -1.25, box[2] - ext)};
  const T hi[2] = {lpt_min((T) 1.25, box[1] + ext),
                   lpt_min((T) 1.25, box[3] + ext)};
  const T xmax = lpt_max(fabs(x0.x), lpt_max(fabs(x0.y), fabs(x0.z)));
  const T tolf = (T) 64 * eps * (rec[LPT_FR_SCALE] + xmax);
  const T tolh = (T) 4 * eps * (rec[LPT_FR_SCALE] + xmax);
  int status = LPT_NEWTON_NOT_CONVERGED;
  s = lpt_clamp(s0, lo[0], hi[0]);
  t = lpt_clamp(t0, lo[1], hi[1]);
  alpha = lpt_clamp(alpha0, (T) -0.5, (T) 1.5);
  iters = 0;

  lpt_v3<T> x, xs, xt, ns, nt;
  if (!lpt_frame(g, f, s, t, x, xs, xt, n, exact, ns, nt))
    return LPT_NEWTON_SINGULAR;
  lpt_v3<T> fv = x0 + alpha * dx - x + a * n;
  T fn = nrm3(fv);

  for (int it = 1; it <= LPT_NEWTON_MAX_ITER; it++) {
    iters = it;
    if (fn <= tolh) return LPT_NEWTON_CONVERGED;

    lpt_v3<T> c1 = dx, c2 = (T) -1 * xs, c3 = (T) -1 * xt;
    if (exact && a > (T) 0) {
      c2 = c2 + a * ns;
      c3 = c3 + a * nt;
    }
    const lpt_v3<T> c23 = cross3(c2, c3);
    const T det = dot3(c1, c23);
    const T dscale = nrm3(c1) * nrm3(c2) * nrm3(c3);
    if (fabs(det) <= (T) 1.0e3 * eps * dscale)
      return fn <= (T) 16 * tolf ? LPT_NEWTON_CONVERGED : LPT_NEWTON_SINGULAR;
    const T d1 = -dot3(fv, c23) / det;
    const T d2 = -dot3(c1, cross3(fv, c3)) / det;
    const T d3 = -dot3(c1, cross3(c2, fv)) / det;

    T lam = 1, atry = alpha, stry = s, ttry = t, fntry = lpt_huge<T>();
    lpt_v3<T> xtry, xstry, xttry, ntry, nstry, nttry, ftry;
    bool oktry = false;
    for (int ih = 1; ih <= LPT_NEWTON_MAX_HALVINGS; ih++) {
      atry = lpt_clamp(alpha + lam * d1, (T) -0.5, (T) 1.5);
      stry = lpt_clamp(s + lam * d2, lo[0], hi[0]);
      ttry = lpt_clamp(t + lam * d3, lo[1], hi[1]);
      oktry = lpt_frame(g, f, stry, ttry, xtry, xstry, xttry, ntry, exact,
                        nstry, nttry);
      if (oktry) {
        ftry = x0 + atry * dx - xtry + a * ntry;
        fntry = nrm3(ftry);
        if (fntry < fn) break;
      }
      lam = (T) 0.5 * lam;
    }
    if (!oktry || fntry >= fn) {
      if (fn <= (T) 16 * tolf) status = LPT_NEWTON_CONVERGED;
      else if (lpt_at_bound(s, t, alpha, lo, hi)) status = LPT_NEWTON_OUTSIDE;
      return status;
    }
    const T stepn = lpt_max(fabs(atry - alpha),
                            lpt_max(fabs(stry - s), fabs(ttry - t)));
    alpha = atry; s = stry; t = ttry;
    x = xtry; xs = xstry; xt = xttry; n = ntry; ns = nstry; nt = nttry;
    fv = ftry; fn = fntry;
    if (stepn <= (T) 100 * eps && fn <= (T) 16 * tolf)
      return LPT_NEWTON_CONVERGED;
  }
  if (fn <= (T) 16 * tolf) status = LPT_NEWTON_CONVERGED;
  else if (lpt_at_bound(s, t, alpha, lo, hi)) status = LPT_NEWTON_OUTSIDE;
  return status;
}


__device__ inline void lpt_edge_st(int k, double u, double &s, double &t) {
  switch (k) {
  case 1: s = -1; t = u; break;
  case 2: s = 1; t = u; break;
  case 3: s = u; t = -1; break;
  default: s = u; t = 1; break;
  }
}

template<typename T>
__device__ inline void lpt_edge_point(const lpt_wall_geom<T> &g, int f, int k, T u,
                            lpt_v3<T> &e, lpt_v3<T> &ep, lpt_v3<T> &epp) {
  double s, t;
  lpt_edge_st(k, (double) u, s, t);
  lpt_v3<T> xs, xt, xss, xst, xtt;
  lpt_eval2(g, f, (T) s, (T) t, e, xs, xt, xss, xst, xtt);
  if (k <= 2) { ep = xt; epp = xtt; }
  else { ep = xs; epp = xss; }
}

template<typename T>
__device__ inline lpt_v3<T> lpt_edge_node(const lpt_wall_geom<T> &g, int f, int k,
                                int i) {
  switch (k) {
  case 1: return xf_node(g, f, 1, i);
  case 2: return xf_node(g, f, g.lx, i);
  case 3: return xf_node(g, f, i, 1);
  default: return xf_node(g, f, i, g.lx);
  }
}

/* segment_segment */
template<typename T>
__device__ inline void lpt_seg_seg(lpt_v3<T> p0, lpt_v3<T> d0, lpt_v3<T> q0,
                         lpt_v3<T> d1, T &s, T &t, T &dist) {
  const T eps = lpt_eps<T>();
  const lpt_v3<T> w = p0 - q0;
  const T a = dot3(d0, d0), b = dot3(d0, d1), c = dot3(d1, d1);
  const T d = dot3(d0, w), e = dot3(d1, w);
  const T den = a * c - b * b;
  s = den > eps * a * c ? lpt_clamp((b * e - c * d) / den, (T) 0, (T) 1)
                        : (T) 0;
  t = c > (T) 0 ? (b * s + e) / c : (T) 0;
  if (t < (T) 0 || t > (T) 1) {
    t = lpt_clamp(t, (T) 0, (T) 1);
    s = a > (T) 0 ? lpt_clamp((b * t - d) / a, (T) 0, (T) 1) : (T) 0;
  }
  dist = nrm3(p0 + s * d0 - q0 - t * d1);
}

/* rim_intersect, x0 facet-local */
template<typename T>
__device__ inline void lpt_rim(const lpt_wall_geom<T> &g, int f, lpt_v3<T> x0,
                     lpt_v3<T> dx, T a, T &alpha, bool &found,
                     lpt_v3<T> &nc, lpt_wall_cnt &cnt) {
  const T eps = lpt_eps<T>();
  const T *rec = fr_of(g, f);
  const T xmax = lpt_max(fabs(x0.x), lpt_max(fabs(x0.y), fabs(x0.z)));
  const T tolf = (T) 64 * eps * (rec[LPT_FR_SCALE] + xmax);

  for (int k = 1; k <= 4; k++) {
    const lpt_v3<T> e0 = fr3(g, f, LPT_FR_CX + 3 * (k - 1));
    const lpt_v3<T> nf = fr3(g, f, LPT_FR_CN + 3 * (k - 1));
    const T qa = dot3(dx, dx);
    const T qb = (T) 2 * dot3(dx, x0 - e0);
    const T qc = dot3(x0 - e0, x0 - e0) - a * a;
    if (qb >= (T) 0 || qa <= (T) 0) continue;
    const T disc = qb * qb - (T) 4 * qa * qc;
    if (disc < (T) 0) continue;
    const T al = lpt_max((T) 0, (-qb - sqrt(disc)) / ((T) 2 * qa));
    if (al > (T) 1 || al >= alpha) continue;
    const lpt_v3<T> x = x0 + al * dx;
    if (dot3(x - e0, nf) > tolf) continue;
    alpha = al;
    found = true;
    nc = ((T) 1 / lpt_max(nrm3(e0 - x), (T) 1e-30)) * (e0 - x);
  }

  for (int k = 1; k <= 4; k++) {
    if (rec[LPT_FR_ELEN + k - 1] <= (T) 0) continue;
    T qa, qb;
    if (!lpt_seg_sphere(x0, dx, fr3(g, f, LPT_FR_EC + 3 * (k - 1)),
                        rec[LPT_FR_ER + k - 1] + a, qa, qb) || qa > alpha)
      continue;

    T dm = lpt_huge<T>(), al_close = 0, u_close = 0;
    lpt_v3<T> ed = mk3<T>(0, 0, 0);
    for (int i = 1; i <= g.lx - 1; i++) {
      const lpt_v3<T> p = lpt_edge_node(g, f, k, i);
      const lpt_v3<T> q = lpt_edge_node(g, f, k, i + 1);
      T al, v, b1;
      lpt_seg_seg(x0, dx, p, q - p, al, v, b1);
      if (b1 < dm) {
        dm = b1;
        al_close = al;
        u_close = g.zg[i - 1] + v * (g.zg[i] - g.zg[i - 1]);
        ed = ((T) 1 / lpt_max(nrm3(q - p), (T) 1e-30)) * (q - p);
      }
    }
    if (dm > a + rec[LPT_FR_EB + k - 1]) continue;
    const lpt_v3<T> dperp = dx - dot3(dx, ed) * ed;
    if (nrm3(dperp) <= eps * nrm3(dx)) continue;
    const T al_entry = al_close -
      sqrt(lpt_max(a * a - dm * dm, (T) 0)) / nrm3(dperp);

    int status = LPT_NEWTON_NOT_CONVERGED;
    bool conv = false;
    T al = 0, u = 0;
    lpt_v3<T> rr;
    for (int attempt = 1; attempt <= 2; attempt++) {
      u = u_close;
      al = lpt_clamp(attempt == 1 ? al_entry : al_close, (T) -0.5, (T) 1.5);
      conv = false;
      status = LPT_NEWTON_NOT_CONVERGED;
      lpt_v3<T> e, ep, epp;
      lpt_edge_point(g, f, k, u, e, ep, epp);
      rr = x0 + al * dx - e;
      T g1 = (dot3(rr, rr) - a * a) / ((T) 2 * a);
      T g2 = dot3(rr, ep) / nrm3(ep);
      T gn = sqrt(g1 * g1 + g2 * g2);
      int it;
      for (it = 1; it <= LPT_NEWTON_MAX_ITER; it++) {
        if (gn <= tolf) { conv = true; break; }
        if (al <= (T) -0.5 || al >= (T) 1.5) {
          status = LPT_NEWTON_OUTSIDE;
          break;
        }
        const T j11 = dot3(rr, dx) / a;
        const T j12 = -dot3(rr, ep) / a;
        const T j21 = dot3(dx, ep) / nrm3(ep);
        const T j22 = (-dot3(ep, ep) + dot3(rr, epp)) / nrm3(ep);
        const T det = j11 * j22 - j12 * j21;
        if (fabs(det) <= (T) 1.0e3 * eps * (fabs(j11 * j22) + fabs(j12 * j21))) {
          /* path tangent to the tube while outside it: no entry */
          if (g1 > (T) 0) status = LPT_NEWTON_OUTSIDE;
          break;
        }
        const T dal = -(j22 * g1 - j12 * g2) / det;
        const T du = -(j11 * g2 - j21 * g1) / det;
        T lam = 1, altry = al, utry = u, g1t = g1, g2t = g2, gnt = gn;
        lpt_v3<T> e1, ept, eppt, rtry = rr;
        for (int ih = 1; ih <= LPT_NEWTON_MAX_HALVINGS; ih++) {
          altry = lpt_clamp(al + lam * dal, (T) -0.5, (T) 1.5);
          utry = lpt_clamp(u + lam * du, (T) -1.25, (T) 1.25);
          lpt_edge_point(g, f, k, utry, e1, ept, eppt);
          rtry = x0 + altry * dx - e1;
          g1t = (dot3(rtry, rtry) - a * a) / ((T) 2 * a);
          g2t = dot3(rtry, ept) / nrm3(ept);
          gnt = sqrt(g1t * g1t + g2t * g2t);
          if (gnt < gn) break;
          lam = (T) 0.5 * lam;
        }
        if (gnt >= gn) {
          conv = gn <= (T) 16 * tolf;
          if (!conv && (g1 > (T) 0 || fabs(u) >= (T) 1.25))
            status = LPT_NEWTON_OUTSIDE;
          break;
        }
        al = altry; u = utry; ep = ept; epp = eppt; rr = rtry;
        g1 = g1t; g2 = g2t; gn = gnt;
      }
      if (conv) status = LPT_NEWTON_CONVERGED;
      if (status != LPT_NEWTON_NOT_CONVERGED) break;
    }
    if (status == LPT_NEWTON_NOT_CONVERGED) cnt.c[LPT_ST_UNRESOLVED] += 1;
    if (!conv) continue;
    if (fabs(u) > (T) 1 || al < (T) 0 || al > (T) 1 || al >= alpha) continue;
    if (dot3(dx, rr) >= (T) 0) continue;
    double s, t;
    lpt_edge_st(k, (double) u, s, t);
    lpt_v3<T> x, xs, xt, nfe, ns, nt;
    if (!lpt_frame(g, f, (T) s, (T) t, x, xs, xt, nfe, false, ns, nt) ||
        dot3(rr, nfe) > tolf)
      continue;
    alpha = al;
    found = true;
    nc = ((T) -1 / nrm3(rr)) * rr;
  }
}

/* Record a converged root, see accept_root in lpt_wall_facet_intersect */
template<typename T>
__device__ inline bool lpt_accept(T s, T t, T al, lpt_v3<T> n, lpt_v3<T> dx,
                        T tol_st, T tol_al, T &alpha, bool &found,
                        lpt_v3<T> &nc) {
  const bool ok = fabs(s) <= (T) 1 + tol_st && fabs(t) <= (T) 1 + tol_st &&
    al >= -tol_al && al <= (T) 1 && dot3(dx, n) > (T) 0;
  if (ok && lpt_max((T) 0, al) < alpha) {
    alpha = lpt_max((T) 0, al);
    found = true;
    nc = n;
  }
  return ok;
}

/* lpt_wall_facet_intersect, x0g in global coordinates */
template<typename T>
__device__ inline void lpt_facet_intersect(const lpt_wall_geom<T> &g, int f,
                                 lpt_v3<T> x0g, lpt_v3<T> dx, T a,
                                 T sg, T tg, bool have_guess, T &alpha,
                                 bool &found, lpt_v3<T> &nc, bool &is_rim,
                                 T tie, lpt_wall_cnt &cnt) {
  const T eps = lpt_eps<T>();
  const T pi = (T) 3.14159265358979323846;
  const T *rec = fr_of(g, f);
  found = false;
  is_rim = false;
  nc = mk3<T>(0, 0, 0);
  const T dxl = nrm3(dx);
  if (!(dxl > (T) 0)) return;
  const lpt_v3<T> dxn = ((T) 1 / dxl) * dx;
  const T tol_st = (T) 1.0e3 * eps;
  const T tol_al = (T) 64 * eps * rec[LPT_FR_SCALE] / dxl;
  const lpt_v3<T> x0 = x0g - fr3(g, f, LPT_FR_ORIGIN);
  const lpt_v3<T> cl = fr3(g, f, LPT_FR_CENTER);
  const lpt_v3<T> ax0 = fr3(g, f, LPT_FR_AXIS);

  T a_in, a_out;
  if (!lpt_seg_sphere(x0, dx, cl, rec[LPT_FR_RADIUS] + a, a_in, a_out) ||
      a_in > alpha)
    return;
  T hmax = lpt_max(dot3(x0 - cl, ax0), dot3(x0 + dx - cl, ax0));
  if (hmax < rec[LPT_FR_SLAB] - a) return;

  T al_r = alpha;
  bool fd_r = false;
  lpt_v3<T> nc_r = mk3<T>(0, 0, 0);
  if (a > (T) 0) lpt_rim(g, f, x0, dx, a, al_r, fd_r, nc_r, cnt);

  T boxes[LPT_PATCH_STACK][4];
  int depth[LPT_PATCH_STACK];
  int top = 0;
  boxes[0][0] = -1; boxes[0][1] = 1; boxes[0][2] = -1; boxes[0][3] = 1;
  depth[0] = 0;
  while (top >= 0) {
    T box[4] = {boxes[top][0], boxes[top][1], boxes[top][2], boxes[top][3]};
    const int dep = depth[top];
    top--;

    lpt_v3<T> c, axis;
    T r, theta, slab[2];
    if (dep == 0) {
      c = cl; r = rec[LPT_FR_RADIUS]; axis = ax0; theta = rec[LPT_FR_ANGLE];
      slab[0] = rec[LPT_FR_SLAB]; slab[1] = rec[LPT_FR_SLAB + 1];
    } else {
      lpt_patch_bounds(g, f, box, c, r, axis, theta, slab);
    }
    if (!lpt_seg_sphere(x0, dx, c, r + a, a_in, a_out) || a_in > alpha)
      continue;
    hmax = lpt_max(dot3(x0 - c, axis), dot3(x0 + dx - c, axis));
    if (hmax < slab[0] - a) continue;

    const T cosang = dot3(dxn, axis);
    bool graph = false;
    if (theta < (T) 0.5 * pi) {
      const T sth = sin(theta);
      if (cosang < -sth) continue;
      graph = cosang > sth;
    }

    bool subdivide = false;
    const T sm = (T) 0.5 * (box[0] + box[1]);
    const T tm = (T) 0.5 * (box[2] + box[3]);
    if (graph || dep == LPT_PATCH_MAX_DEPTH) {
      T s0 = sm, t0 = tm;
      if (dep == 0 && have_guess) {
        s0 = lpt_clamp(sg, (T) -1, (T) 1);
        t0 = lpt_clamp(tg, (T) -1, (T) 1);
      }
      lpt_v3<T> xg, xs, xt, ng, ns, nt;
      T al0 = (T) 0.5 * (a_in + a_out);
      if (lpt_frame(g, f, s0, t0, xg, xs, xt, ng, false, ns, nt)) {
        const T den = dot3(dx, ng);
        if (fabs(den) > eps * dxl)
          al0 = dot3(xg - a * ng - x0, ng) / den;
      }
      T s, t, al;
      lpt_v3<T> nrm;
      int iters;
      int status = lpt_newton(g, f, x0, dx, a, box, s0, t0, al0, false, s,
                              t, al, nrm, iters);
      if ((status == LPT_NEWTON_SINGULAR ||
           status == LPT_NEWTON_NOT_CONVERGED) && a > (T) 0 &&
          rec[LPT_FR_PLANAR] == (T) 0) {
        status = lpt_newton(g, f, x0, dx, a, box, s0, t0, al0, true, s, t,
                            al, nrm, iters);
      }
      bool valid = false;
      if (status == LPT_NEWTON_CONVERGED)
        valid = lpt_accept(s, t, al, nrm, dx, tol_st, tol_al, alpha, found,
                           nc);
      if (!valid && !graph && dep == LPT_PATCH_MAX_DEPTH) {
        status = lpt_newton(g, f, x0, dx, a, box, sm, tm, a_out, a > (T) 0,
                            s, t, al, nrm, iters);
        if (status == LPT_NEWTON_CONVERGED)
          valid = lpt_accept(s, t, al, nrm, dx, tol_st, tol_al, alpha,
                             found, nc);
        if (!valid) cnt.c[LPT_ST_UNRESOLVED] += 1;
      } else if (!valid) {
        if (status == LPT_NEWTON_CONVERGED ||
            (status == LPT_NEWTON_OUTSIDE && graph)) {
          /* the unique crossing is outside this patch or this step */
        } else if (dep < LPT_PATCH_MAX_DEPTH) {
          subdivide = true;
        } else {
          cnt.c[LPT_ST_UNRESOLVED] += 1;
        }
      }
    } else {
      subdivide = true;
    }

    if (subdivide) {
      const T q[4][4] = {{box[0], sm, box[2], tm}, {sm, box[1], box[2], tm},
                         {box[0], sm, tm, box[3]}, {sm, box[1], tm, box[3]}};
      for (int k = 0; k < 4; k++) {
        top++;
        depth[top] = dep + 1;
        for (int i = 0; i < 4; i++) boxes[top][i] = q[k][i];
      }
    }
  }

  if (fd_r && (!found || al_r < alpha - tie)) {
    alpha = al_r;
    nc = nc_r;
    found = true;
    is_rim = true;
  }
}

/* facet_st_from_rst */
template<typename T>
__device__ inline void lpt_st_from_rst(int facet, T r, T s, T t, T &sf, T &tf) {
  switch (facet) {
  case 1: case 2: sf = s; tf = t; break;
  case 3: case 4: sf = r; tf = t; break;
  default: sf = r; tf = s; break;
  }
  sf = lpt_clamp(sf, (T) -1, (T) 1);
  tf = lpt_clamp(tf, (T) -1, (T) 1);
}

/* lpt_wall_collide: el is 1-based; rst may be null */
template<typename T>
__device__ inline int lpt_collide(const lpt_wall_geom<T> &g, int el, lpt_v3<T> p0,
                        lpt_v3<T> &p1, T a, T e_n,
                        lpt_v3<T> normals[LPT_WALL_MAX_BOUNCES],
                        const T *rst, lpt_wall_cnt &cnt) {
  const T eps = lpt_eps<T>();
  int n_hits = 0;
  if (el < 1 || el > g.nelv) return 0;
  if (g.el_ptr[el - 1] >= g.el_ptr[el]) return 0;
  if (nrm3(p1 - p0) + a > g.el_radius[el - 1]) cnt.c[LPT_ST_REACH] += 1;

  lpt_v3<T> q0 = p0, q1 = p1;
  for (int bounce = 1; bounce <= LPT_WALL_MAX_BOUNCES + 1; bounce++) {
    const lpt_v3<T> dx = q1 - q0;
    if (!(nrm3(dx) > (T) 0)) break;
    T alpha = lpt_huge<T>();
    int fbest = 0;
    bool rim_best = false;
    lpt_v3<T> nbest = mk3<T>(0, 0, 0);
    const T qmax = lpt_max(fabs(q0.x), lpt_max(fabs(q0.y), fabs(q0.z)));
    const T a_tie = ((T) 128 * eps * (qmax + nrm3(dx)) + (T) 1.0e-3 * a) /
      nrm3(dx);
    for (int k = g.el_ptr[el - 1]; k < g.el_ptr[el]; k++) {
      const int f = g.el_facets[k - 1];
      bool have_guess = false;
      T sg = 0, tg = 0;
      if (bounce == 1 && rst != 0 && g.fi[2 * (f - 1)] == el &&
          g.fi[2 * (f - 1) + 1] > 0) {
        lpt_st_from_rst(g.fi[2 * (f - 1) + 1], rst[0], rst[1], rst[2], sg,
                        tg);
        have_guess = true;
      }
      T al_c = rim_best ? alpha + a_tie : alpha;
      bool found, rim;
      lpt_v3<T> nc;
      lpt_facet_intersect(g, f, q0, dx, a, sg, tg, have_guess, al_c, found,
                          nc, rim, a_tie, cnt);
      if (found) {
        if (rim && !rim_best) found = al_c < alpha - a_tie;
        else if (!rim && rim_best) found = al_c < alpha + a_tie;
        else found = al_c < alpha;
      }
      if (found) {
        alpha = al_c;
        fbest = f;
        nbest = nc;
        rim_best = rim;
      }
    }
    if (fbest == 0) break;

    const T tol_off = (T) 16 * eps * fr_of(g, fbest)[LPT_FR_SCALE];
    if (n_hits == LPT_WALL_MAX_BOUNCES) {
      q1 = q0 + alpha * dx - tol_off * nbest;
      break;
    }
    lpt_v3<T> r = ((T) 1 - alpha) * dx;
    q0 = q0 + alpha * dx - tol_off * nbest;
    const T rn = dot3(r, nbest);
    if (rn > (T) 0) r = r - ((T) 1 + e_n) * rn * nbest;
    q1 = q0 + r;
    normals[n_hits] = nbest;
    n_hits++;
  }
  p1 = q1;
  return n_hits;
}

/* v <- v - (1 + e_n)(v . n) n */
template<typename T>
__device__ inline void lpt_reflect(T &vx, T &vy, T &vz, lpt_v3<T> n, T e_n) {
  const T vn = ((T) 1 + e_n) * (vx * n.x + vy * n.y + vz * n.z);
  vx -= vn * n.x;
  vy -= vn * n.y;
  vz -= vn * n.z;
}

/* Body of the per-particle kernel, shared by CUDA, HIP and host tests */
template<typename T>
__device__ inline void lpt_wall_collide_particle(
    const lpt_wall_geom<T> &g, int p, const int *el_list, const T *rst,
    const T *x_old, const T *y_old, const T *z_old, T *x, T *y, T *z,
    const T *d, T *u, T *v, T *w, T *u_lag, T *v_lag, T *w_lag,
    T *u_laglag, T *v_laglag, T *w_laglag, T *acc_xlag, T *acc_ylag,
    T *acc_zlag, T *acc_xlaglag, T *acc_ylaglag, T *acc_zlaglag,
    T *u_old, T *v_old, T *w_old, T *acc_x, T *acc_y, T *acc_z,
    int lag_len, T e_n, lpt_wall_cnt &cnt) {
  if (el_list[p] < 0) return;
  lpt_v3<T> p1 = mk3(x[p], y[p], z[p]);
  lpt_v3<T> normals[LPT_WALL_MAX_BOUNCES];
  const int n_hits = lpt_collide(g, el_list[p] + 1,
                                 mk3(x_old[p], y_old[p], z_old[p]), p1,
                                 (T) 0.5 * d[p], e_n, normals,
                                 rst != 0 ? rst + 3 * (size_t) p : (const T *) 0,
                                 cnt);
  x[p] = p1.x;
  y[p] = p1.y;
  z[p] = p1.z;
  for (int k = 0; k < n_hits; k++) {
    const lpt_v3<T> nn = normals[k];
    lpt_reflect(u[p], v[p], w[p], nn, e_n);
    lpt_reflect(u_old[p], v_old[p], w_old[p], nn, e_n);
    lpt_reflect(acc_x[p], acc_y[p], acc_z[p], nn, e_n);
    if (lag_len >= 1) {
      lpt_reflect(u_lag[p], v_lag[p], w_lag[p], nn, e_n);
      lpt_reflect(acc_xlag[p], acc_ylag[p], acc_zlag[p], nn, e_n);
    }
    if (lag_len >= 2) {
      lpt_reflect(u_laglag[p], v_laglag[p], w_laglag[p], nn, e_n);
      lpt_reflect(acc_xlaglag[p], acc_ylaglag[p], acc_zlaglag[p], nn, e_n);
    }
  }
}

#endif /* __LPT_WALL_INTERSECT_H__ */
