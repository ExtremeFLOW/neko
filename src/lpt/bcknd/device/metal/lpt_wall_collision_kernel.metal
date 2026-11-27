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
 * Metal kernels for elastic collisions of particles with the high-order
 * wall facets, see lpt_wall_geometry.f90 for the formulation. The shared
 * part is identical to the OpenCL kernel (with float literals).
 *
 * Metal allows at most 31 buffers per function, so the work is split:
 * lpt_wall_collide_kernel finds the collisions, corrects the positions and
 * stores the contact normals of each particle; lpt_wall_reflect_kernel then
 * applies them to the velocity and its history.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

typedef float real;
typedef float3 real3;
static inline real3 lw_v3(real a, real b, real c) { return real3(a, b, c); }

/*
 * Trajectory-based collision of particles with the high-order wall facets.
 * C-style port of lpt_wall_geometry.f90 and of the CUDA/HIP header
 * lpt/bcknd/device/lpt_wall_intersect.h. The OpenCL kernel holds a separate
 * copy of this code; keep the two in sync.
 *
 * The per-facet record layout and the statistics counters are those of
 * lpt_wall_intersect.h. Indices of facets and elements are 1-based.
 */

#define LPT_WALL_MAX_LX 16
#define LPT_WALL_MAX_BOUNCES 8
#define LPT_NEWTON_MAX_ITER 25
#define LPT_NEWTON_MAX_HALVINGS 10
#define LPT_PATCH_MAX_DEPTH 4
#define LPT_PATCH_STACK (3 * LPT_PATCH_MAX_DEPTH + 2)

#define LPT_NEWTON_CONVERGED 0
#define LPT_NEWTON_SINGULAR 1
#define LPT_NEWTON_NOT_CONVERGED 2
#define LPT_NEWTON_OUTSIDE 3

#define LPT_FR_ORIGIN 0
#define LPT_FR_ORIENT 3
#define LPT_FR_CENTER 4
#define LPT_FR_RADIUS 7
#define LPT_FR_AXIS 8
#define LPT_FR_ANGLE 11
#define LPT_FR_SLAB 12
#define LPT_FR_SCALE 14
#define LPT_FR_PLANAR 15
#define LPT_FR_ELEN 16
#define LPT_FR_EC 20
#define LPT_FR_ER 32
#define LPT_FR_EB 36
#define LPT_FR_CX 40
#define LPT_FR_CN 52
#define LPT_NFR 64

#define LPT_ST_UNRESOLVED 0
#define LPT_ST_REACH 1
#define LPT_NSTATS 2

/* Machine epsilon of the working precision, 2^-23 or 2^-52 */
#define LPT_EPS (sizeof(real) == 8 ? (real) 2.220446049250313e-16f \
                 : (real) 1.1920928955078125e-7f)
#define LPT_HUGE ((real) 1.0e30f)

typedef struct {
  device const real *xf;
  device const real *fr;
  device const int *fi;
  device const real *zg;
  device const int *el_ptr;
  device const int *el_facets;
  device const real *el_radius;
  int lx;
  int nelv;
} lw_geom;

static inline real lw_max(real a, real b) { return a > b ? a : b; }
static inline real lw_min(real a, real b) { return a < b ? a : b; }
static inline real lw_clamp(real v, real lo, real hi) {
  return lw_max(lo, lw_min(hi, v));
}
static inline real lw_dot(real3 a, real3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline real lw_nrm(real3 a) { return sqrt(lw_dot(a, a)); }
static inline real3 lw_cross(real3 a, real3 b) {
  return lw_v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
}
static inline real3 lw_scale(real s, real3 a) {
  return lw_v3(s * a.x, s * a.y, s * a.z);
}
static inline real3 lw_add(real3 a, real3 b) {
  return lw_v3(a.x + b.x, a.y + b.y, a.z + b.z);
}
static inline real3 lw_sub(real3 a, real3 b) {
  return lw_v3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline real lw_fr(lw_geom g, int f, int off) {
  return g.fr[LPT_NFR * (f - 1) + off];
}
static inline real3 lw_fr3(lw_geom g, int f, int off) {
  const int i = LPT_NFR * (f - 1) + off;
  return lw_v3(g.fr[i], g.fr[i + 1], g.fr[i + 2]);
}

/* fd_weights_full; c[3 * j + k], j node, k derivative order */
static inline void lw_fd_weights(real xi, device const real *zg, int lx, int m,
                              thread real *c) {
  const int n = lx - 1;
  real c1 = 1.0f;
  real c4 = zg[0] - xi;
  for (int j = 0; j <= n; j++)
    for (int k = 0; k <= m; k++)
      c[3 * j + k] = 0.0f;
  c[0] = 1.0f;
  for (int i = 1; i <= n; i++) {
    const int mn = i < m ? i : m;
    real c2 = 1.0f;
    const real c5 = c4;
    c4 = zg[i] - xi;
    for (int j = 0; j <= i - 1; j++) {
      const real c3 = zg[i] - zg[j];
      c2 = c2 * c3;
      for (int k = mn; k >= 1; k--)
        c[3 * i + k] = c1 * (k * c[3 * (i - 1) + k - 1]
                             - c5 * c[3 * (i - 1) + k]) / c2;
      c[3 * i] = -c1 * c5 * c[3 * (i - 1)] / c2;
      for (int k = mn; k >= 1; k--)
        c[3 * j + k] = (c4 * c[3 * j + k] - k * c[3 * j + k - 1]) / c3;
      c[3 * j] = c4 * c[3 * j] / c3;
    }
    c1 = c2;
  }
}

static inline real3 lw_node(lw_geom g, int f, int a, int b) {
  const int i = 3 * ((a - 1) + g.lx * ((b - 1) + g.lx * (f - 1)));
  return lw_v3(g.xf[i], g.xf[i + 1], g.xf[i + 2]);
}

/* eval_local / eval_nodes */
static inline void lw_eval2(lw_geom g, int f, real s, real t, thread real3 *x,
                         thread real3 *xs, thread real3 *xt, thread real3 *xss,
                         thread real3 *xst, thread real3 *xtt) {
  real cs[3 * LPT_WALL_MAX_LX], ct[3 * LPT_WALL_MAX_LX];
  lw_fd_weights(s, g.zg, g.lx, 2, cs);
  lw_fd_weights(t, g.zg, g.lx, 2, ct);
  const real3 z = lw_v3(0.0f, 0.0f, 0.0f);
  *x = z; *xs = z; *xt = z; *xss = z; *xst = z; *xtt = z;
  for (int b = 1; b <= g.lx; b++) {
    real3 a0 = z, a1 = z, a2 = z;
    for (int a = 1; a <= g.lx; a++) {
      const real3 p = lw_node(g, f, a, b);
      a0 = lw_add(a0, lw_scale(cs[3 * (a - 1)], p));
      a1 = lw_add(a1, lw_scale(cs[3 * (a - 1) + 1], p));
      a2 = lw_add(a2, lw_scale(cs[3 * (a - 1) + 2], p));
    }
    *x = lw_add(*x, lw_scale(ct[3 * (b - 1)], a0));
    *xs = lw_add(*xs, lw_scale(ct[3 * (b - 1)], a1));
    *xt = lw_add(*xt, lw_scale(ct[3 * (b - 1) + 1], a0));
    *xss = lw_add(*xss, lw_scale(ct[3 * (b - 1)], a2));
    *xst = lw_add(*xst, lw_scale(ct[3 * (b - 1) + 1], a1));
    *xtt = lw_add(*xtt, lw_scale(ct[3 * (b - 1) + 2], a0));
  }
}

/* eval_first */
static inline void lw_eval1(lw_geom g, int f, real s, real t, thread real3 *x,
                         thread real3 *xs, thread real3 *xt) {
  real cs[3 * LPT_WALL_MAX_LX], ct[3 * LPT_WALL_MAX_LX];
  lw_fd_weights(s, g.zg, g.lx, 1, cs);
  lw_fd_weights(t, g.zg, g.lx, 1, ct);
  const real3 z = lw_v3(0.0f, 0.0f, 0.0f);
  *x = z; *xs = z; *xt = z;
  for (int b = 1; b <= g.lx; b++) {
    real3 a0 = z, a1 = z;
    for (int a = 1; a <= g.lx; a++) {
      const real3 p = lw_node(g, f, a, b);
      a0 = lw_add(a0, lw_scale(cs[3 * (a - 1)], p));
      a1 = lw_add(a1, lw_scale(cs[3 * (a - 1) + 1], p));
    }
    *x = lw_add(*x, lw_scale(ct[3 * (b - 1)], a0));
    *xs = lw_add(*xs, lw_scale(ct[3 * (b - 1)], a1));
    *xt = lw_add(*xt, lw_scale(ct[3 * (b - 1) + 1], a0));
  }
}

/* surface_frame; returns 0 if the normal is degenerate */
static inline int lw_frame(lw_geom g, int f, real s, real t, thread real3 *x,
                        thread real3 *xs, thread real3 *xt, thread real3 *n,
                        int deriv, thread real3 *ns, thread real3 *nt) {
  real3 xss, xst, xtt;
  if (deriv)
    lw_eval2(g, f, s, t, x, xs, xt, &xss, &xst, &xtt);
  else
    lw_eval1(g, f, s, t, x, xs, xt);
  const real3 nv = lw_cross(*xs, *xt);
  const real nlen = lw_nrm(nv);
  const real orient = lw_fr(g, f, LPT_FR_ORIENT);
  if (!(nlen > 0.0f)) {
    *n = lw_v3(0.0f, 0.0f, 0.0f); *ns = *n; *nt = *n;
    return 0;
  }
  const real3 nu = lw_scale(1.0f / nlen, nv);
  *n = lw_scale(orient, nu);
  if (deriv) {
    const real3 ds = lw_add(lw_cross(xss, *xt), lw_cross(*xs, xst));
    const real3 dt = lw_add(lw_cross(xst, *xt), lw_cross(*xs, xtt));
    *ns = lw_scale(orient / nlen, lw_sub(ds, lw_scale(lw_dot(nu, ds), nu)));
    *nt = lw_scale(orient / nlen, lw_sub(dt, lw_scale(lw_dot(nu, dt), nu)));
  } else {
    *ns = lw_v3(0.0f, 0.0f, 0.0f); *nt = *ns;
  }
  return 1;
}

static inline real lw_angle(real3 u, real3 v) {
  return atan2(lw_nrm(lw_cross(u, v)), lw_dot(u, v));
}

/* patch_bounds, two passes over the samples */
static inline void lw_patch_bounds(lw_geom g, int f, thread const real *box,
                                thread real3 *c, thread real *r,
                                thread real3 *axis, thread real *theta,
                                thread real *slab) {
  const int m = 2 * g.lx - 1 > 3 ? 2 * g.lx - 1 : 3;
  const real hs = (box[1] - box[0]) / (m - 1);
  const real ht = (box[3] - box[2]) / (m - 1);
  const real orient = lw_fr(g, f, LPT_FR_ORIENT);
  real3 sx = lw_v3(0.0f, 0.0f, 0.0f), sn = sx, x, xs, xt;
  real gs = 0.0f, gt = 0.0f;
  int degenerate = 0;
  for (int b = 1; b <= m; b++) {
    for (int a = 1; a <= m; a++) {
      lw_eval1(g, f, box[0] + (a - 1) * hs, box[2] + (b - 1) * ht, &x, &xs,
               &xt);
      gs = lw_max(gs, lw_nrm(xs));
      gt = lw_max(gt, lw_nrm(xt));
      const real3 nv = lw_cross(xs, xt);
      const real nl = lw_nrm(nv);
      if (nl > 0.0f) sn = lw_add(sn, lw_scale(orient / nl, nv));
      else degenerate = 1;
      sx = lw_add(sx, x);
    }
  }
  *c = lw_scale(1.0f / (real) (m * m), sx);
  const real margin = 0.5f * (hs * gs + ht * gt);
  const real nl = lw_nrm(sn);
  *r = 0.0f;
  if (degenerate || !(nl > 0.0f)) {
    for (int b = 1; b <= m; b++)
      for (int a = 1; a <= m; a++) {
        lw_eval1(g, f, box[0] + (a - 1) * hs, box[2] + (b - 1) * ht, &x,
                 &xs, &xt);
        *r = lw_max(*r, lw_nrm(lw_sub(x, *c)));
      }
    *r = 1.01f * (*r + margin);
    *axis = lw_v3(0.0f, 0.0f, 1.0f);
    *theta = LPT_HUGE;
    slab[0] = -*r;
    slab[1] = *r;
    return;
  }
  *axis = lw_scale(1.0f / nl, sn);
  *theta = 0.0f;
  real dmax = 0.0f;
  slab[0] = LPT_HUGE;
  slab[1] = -LPT_HUGE;
  real3 prev_row[2 * LPT_WALL_MAX_LX + 1];
  for (int b = 1; b <= m; b++) {
    real3 left = lw_v3(0.0f, 0.0f, 0.0f);
    for (int a = 1; a <= m; a++) {
      lw_eval1(g, f, box[0] + (a - 1) * hs, box[2] + (b - 1) * ht, &x, &xs,
               &xt);
      *r = lw_max(*r, lw_nrm(lw_sub(x, *c)));
      const real h = lw_dot(lw_sub(x, *c), *axis);
      slab[0] = lw_min(slab[0], h);
      slab[1] = lw_max(slab[1], h);
      const real3 nv = lw_cross(xs, xt);
      const real3 nn = lw_scale(orient / lw_nrm(nv), nv);
      *theta = lw_max(*theta, lw_angle(nn, *axis));
      if (a > 1) dmax = lw_max(dmax, lw_angle(nn, left));
      if (b > 1) dmax = lw_max(dmax, lw_angle(nn, prev_row[a - 1]));
      left = nn;
      prev_row[a - 1] = nn;
    }
  }
  *r = 1.01f * (*r + margin);
  *theta = 1.01f * (*theta + dmax);
  slab[0] -= 1.01f * margin;
  slab[1] += 1.01f * margin;
}

/* segment_sphere_interval */
static inline int lw_seg_sphere(real3 x0, real3 dx, real3 c, real r,
                             thread real *a_in, thread real *a_out) {
  *a_in = 0.0f;
  *a_out = 0.0f;
  const real3 w = lw_sub(x0, c);
  const real qa = lw_dot(dx, dx);
  const real qb = lw_dot(dx, w);
  const real qc = lw_dot(w, w) - r * r;
  if (qa <= 0.0f) {
    *a_out = 1.0f;
    return qc <= 0.0f;
  }
  const real disc = qb * qb - qa * qc;
  if (disc < 0.0f) return 0;
  const real sq = sqrt(disc);
  *a_in = lw_max(0.0f, (-qb - sq) / qa);
  *a_out = lw_min(1.0f, (-qb + sq) / qa);
  return *a_in <= *a_out;
}

static inline int lw_at_bound(real s, real t, real al, thread const real *lo,
                           thread const real *hi) {
  return s <= lo[0] || t <= lo[1] || s >= hi[0] || t >= hi[1] ||
    al <= -0.5f || al >= 1.5f;
}

/* newton_intersect */
static inline int lw_newton(lw_geom g, int f, real3 x0, real3 dx, real a,
                         thread const real *box, real s0, real t0,
                         real alpha0, int exact, thread real *s,
                         thread real *t, thread real *alpha, thread real3 *n,
                         thread int *iters) {
  const real eps = LPT_EPS;
  const real ext = 0.5f * lw_max(box[1] - box[0], box[3] - box[2]);
  real lo[2], hi[2];
  lo[0] = lw_max(-1.25f, box[0] - ext);
  lo[1] = lw_max(-1.25f, box[2] - ext);
  hi[0] = lw_min(1.25f, box[1] + ext);
  hi[1] = lw_min(1.25f, box[3] + ext);
  const real xmax = lw_max(fabs(x0.x), lw_max(fabs(x0.y), fabs(x0.z)));
  const real tolf = 64.0f * eps * (lw_fr(g, f, LPT_FR_SCALE) + xmax);
  const real tolh = 4.0f * eps * (lw_fr(g, f, LPT_FR_SCALE) + xmax);
  int status = LPT_NEWTON_NOT_CONVERGED;
  *s = lw_clamp(s0, lo[0], hi[0]);
  *t = lw_clamp(t0, lo[1], hi[1]);
  *alpha = lw_clamp(alpha0, -0.5f, 1.5f);
  *iters = 0;

  real3 x, xs, xt, ns, nt;
  if (!lw_frame(g, f, *s, *t, &x, &xs, &xt, n, exact, &ns, &nt))
    return LPT_NEWTON_SINGULAR;
  real3 fv = lw_add(lw_sub(lw_add(x0, lw_scale(*alpha, dx)), x),
                    lw_scale(a, *n));
  real fn = lw_nrm(fv);

  for (int it = 1; it <= LPT_NEWTON_MAX_ITER; it++) {
    *iters = it;
    if (fn <= tolh) return LPT_NEWTON_CONVERGED;

    const real3 c1 = dx;
    real3 c2 = lw_scale(-1.0f, xs), c3 = lw_scale(-1.0f, xt);
    if (exact && a > 0.0f) {
      c2 = lw_add(c2, lw_scale(a, ns));
      c3 = lw_add(c3, lw_scale(a, nt));
    }
    const real3 c23 = lw_cross(c2, c3);
    const real det = lw_dot(c1, c23);
    const real dscale = lw_nrm(c1) * lw_nrm(c2) * lw_nrm(c3);
    if (fabs(det) <= 1.0e3f * eps * dscale)
      return fn <= 16.0f * tolf ? LPT_NEWTON_CONVERGED : LPT_NEWTON_SINGULAR;
    const real d1 = -lw_dot(fv, c23) / det;
    const real d2 = -lw_dot(c1, lw_cross(fv, c3)) / det;
    const real d3 = -lw_dot(c1, lw_cross(c2, fv)) / det;

    real lam = 1.0f, atry = *alpha, stry = *s, ttry = *t, fntry = LPT_HUGE;
    real3 xtry, xstry, xttry, ntry, nstry, nttry, ftry;
    int oktry = 0;
    for (int ih = 1; ih <= LPT_NEWTON_MAX_HALVINGS; ih++) {
      atry = lw_clamp(*alpha + lam * d1, -0.5f, 1.5f);
      stry = lw_clamp(*s + lam * d2, lo[0], hi[0]);
      ttry = lw_clamp(*t + lam * d3, lo[1], hi[1]);
      oktry = lw_frame(g, f, stry, ttry, &xtry, &xstry, &xttry, &ntry,
                       exact, &nstry, &nttry);
      if (oktry) {
        ftry = lw_add(lw_sub(lw_add(x0, lw_scale(atry, dx)), xtry),
                      lw_scale(a, ntry));
        fntry = lw_nrm(ftry);
        if (fntry < fn) break;
      }
      lam = 0.5f * lam;
    }
    if (!oktry || fntry >= fn) {
      if (fn <= 16.0f * tolf) status = LPT_NEWTON_CONVERGED;
      else if (lw_at_bound(*s, *t, *alpha, lo, hi))
        status = LPT_NEWTON_OUTSIDE;
      return status;
    }
    const real stepn = lw_max(fabs(atry - *alpha),
                              lw_max(fabs(stry - *s), fabs(ttry - *t)));
    *alpha = atry; *s = stry; *t = ttry;
    x = xtry; xs = xstry; xt = xttry; *n = ntry; ns = nstry; nt = nttry;
    fv = ftry; fn = fntry;
    if (stepn <= 100.0f * eps && fn <= 16.0f * tolf)
      return LPT_NEWTON_CONVERGED;
  }
  if (fn <= 16.0f * tolf) status = LPT_NEWTON_CONVERGED;
  else if (lw_at_bound(*s, *t, *alpha, lo, hi)) status = LPT_NEWTON_OUTSIDE;
  return status;
}


static inline void lw_edge_st(int k, real u, thread real *s, thread real *t) {
  if (k == 1) { *s = -1.0f; *t = u; }
  else if (k == 2) { *s = 1.0f; *t = u; }
  else if (k == 3) { *s = u; *t = -1.0f; }
  else { *s = u; *t = 1.0f; }
}

static inline void lw_edge_point(lw_geom g, int f, int k, real u,
                              thread real3 *e, thread real3 *ep,
                              thread real3 *epp) {
  real s, t;
  lw_edge_st(k, u, &s, &t);
  real3 xs, xt, xss, xst, xtt;
  lw_eval2(g, f, s, t, e, &xs, &xt, &xss, &xst, &xtt);
  if (k <= 2) { *ep = xt; *epp = xtt; }
  else { *ep = xs; *epp = xss; }
}

static inline real3 lw_edge_node(lw_geom g, int f, int k, int i) {
  if (k == 1) return lw_node(g, f, 1, i);
  if (k == 2) return lw_node(g, f, g.lx, i);
  if (k == 3) return lw_node(g, f, i, 1);
  return lw_node(g, f, i, g.lx);
}

/* segment_segment */
static inline void lw_seg_seg(real3 p0, real3 d0, real3 q0, real3 d1,
                           thread real *s, thread real *t, thread real *dist) {
  const real eps = LPT_EPS;
  const real3 w = lw_sub(p0, q0);
  const real a = lw_dot(d0, d0), b = lw_dot(d0, d1), c = lw_dot(d1, d1);
  const real d = lw_dot(d0, w), e = lw_dot(d1, w);
  const real den = a * c - b * b;
  *s = den > eps * a * c ? lw_clamp((b * e - c * d) / den, 0.0f, 1.0f) : 0.0f;
  *t = c > 0.0f ? (b * *s + e) / c : 0.0f;
  if (*t < 0.0f || *t > 1.0f) {
    *t = lw_clamp(*t, 0.0f, 1.0f);
    *s = a > 0.0f ? lw_clamp((b * *t - d) / a, 0.0f, 1.0f) : 0.0f;
  }
  *dist = lw_nrm(lw_sub(lw_add(p0, lw_scale(*s, d0)),
                        lw_add(q0, lw_scale(*t, d1))));
}

/* rim_intersect, x0 facet-local */
static inline void lw_rim(lw_geom g, int f, real3 x0, real3 dx, real a,
                       thread real *alpha, thread int *found,
                       thread real3 *nc, thread uint *cnt) {
  const real eps = LPT_EPS;
  const real xmax = lw_max(fabs(x0.x), lw_max(fabs(x0.y), fabs(x0.z)));
  const real tolf = 64.0f * eps * (lw_fr(g, f, LPT_FR_SCALE) + xmax);

  for (int k = 1; k <= 4; k++) {
    const real3 e0 = lw_fr3(g, f, LPT_FR_CX + 3 * (k - 1));
    const real3 nf = lw_fr3(g, f, LPT_FR_CN + 3 * (k - 1));
    const real3 w = lw_sub(x0, e0);
    const real qa = lw_dot(dx, dx);
    const real qb = 2.0f * lw_dot(dx, w);
    const real qc = lw_dot(w, w) - a * a;
    if (qb >= 0.0f || qa <= 0.0f) continue;
    const real disc = qb * qb - 4.0f * qa * qc;
    if (disc < 0.0f) continue;
    const real al = lw_max(0.0f, (-qb - sqrt(disc)) / (2.0f * qa));
    if (al > 1.0f || al >= *alpha) continue;
    const real3 x = lw_add(x0, lw_scale(al, dx));
    if (lw_dot(lw_sub(x, e0), nf) > tolf) continue;
    *alpha = al;
    *found = 1;
    *nc = lw_scale(1.0f / lw_max(lw_nrm(lw_sub(e0, x)), 1.0e-30f),
                   lw_sub(e0, x));
  }

  for (int k = 1; k <= 4; k++) {
    if (lw_fr(g, f, LPT_FR_ELEN + k - 1) <= 0.0f) continue;
    real qa, qb;
    if (!lw_seg_sphere(x0, dx, lw_fr3(g, f, LPT_FR_EC + 3 * (k - 1)),
                       lw_fr(g, f, LPT_FR_ER + k - 1) + a, &qa, &qb) ||
        qa > *alpha)
      continue;

    real dm = LPT_HUGE, al_close = 0.0f, u_close = 0.0f;
    real3 ed = lw_v3(0.0f, 0.0f, 0.0f);
    for (int i = 1; i <= g.lx - 1; i++) {
      const real3 p = lw_edge_node(g, f, k, i);
      const real3 q = lw_edge_node(g, f, k, i + 1);
      real al, v, b1;
      lw_seg_seg(x0, dx, p, lw_sub(q, p), &al, &v, &b1);
      if (b1 < dm) {
        dm = b1;
        al_close = al;
        u_close = g.zg[i - 1] + v * (g.zg[i] - g.zg[i - 1]);
        ed = lw_scale(1.0f / lw_max(lw_nrm(lw_sub(q, p)), 1.0e-30f),
                      lw_sub(q, p));
      }
    }
    if (dm > a + lw_fr(g, f, LPT_FR_EB + k - 1)) continue;
    const real3 dperp = lw_sub(dx, lw_scale(lw_dot(dx, ed), ed));
    if (lw_nrm(dperp) <= eps * lw_nrm(dx)) continue;
    const real al_entry = al_close -
      sqrt(lw_max(a * a - dm * dm, 0.0f)) / lw_nrm(dperp);

    int status = LPT_NEWTON_NOT_CONVERGED;
    int conv = 0;
    real al = 0.0f, u = 0.0f;
    real3 rr = lw_v3(0.0f, 0.0f, 0.0f);
    for (int attempt = 1; attempt <= 2; attempt++) {
      u = u_close;
      al = lw_clamp(attempt == 1 ? al_entry : al_close, -0.5f, 1.5f);
      conv = 0;
      status = LPT_NEWTON_NOT_CONVERGED;
      real3 e, ep, epp;
      lw_edge_point(g, f, k, u, &e, &ep, &epp);
      rr = lw_sub(lw_add(x0, lw_scale(al, dx)), e);
      real g1 = (lw_dot(rr, rr) - a * a) / (2.0f * a);
      real g2 = lw_dot(rr, ep) / lw_nrm(ep);
      real gn = sqrt(g1 * g1 + g2 * g2);
      int it;
      for (it = 1; it <= LPT_NEWTON_MAX_ITER; it++) {
        if (gn <= tolf) { conv = 1; break; }
        if (al <= -0.5f || al >= 1.5f) {
          status = LPT_NEWTON_OUTSIDE;
          break;
        }
        const real j11 = lw_dot(rr, dx) / a;
        const real j12 = -lw_dot(rr, ep) / a;
        const real j21 = lw_dot(dx, ep) / lw_nrm(ep);
        const real j22 = (-lw_dot(ep, ep) + lw_dot(rr, epp)) / lw_nrm(ep);
        const real det = j11 * j22 - j12 * j21;
        if (fabs(det) <= 1.0e3f * eps * (fabs(j11 * j22) + fabs(j12 * j21))) {
          /* path tangent to the tube while outside it: no entry */
          if (g1 > 0.0f) status = LPT_NEWTON_OUTSIDE;
          break;
        }
        const real dal = -(j22 * g1 - j12 * g2) / det;
        const real du = -(j11 * g2 - j21 * g1) / det;
        real lam = 1.0f, altry = al, utry = u, g1t = g1, g2t = g2, gnt = gn;
        real3 e1, ept, eppt, rtry = rr;
        for (int ih = 1; ih <= LPT_NEWTON_MAX_HALVINGS; ih++) {
          altry = lw_clamp(al + lam * dal, -0.5f, 1.5f);
          utry = lw_clamp(u + lam * du, -1.25f, 1.25f);
          lw_edge_point(g, f, k, utry, &e1, &ept, &eppt);
          rtry = lw_sub(lw_add(x0, lw_scale(altry, dx)), e1);
          g1t = (lw_dot(rtry, rtry) - a * a) / (2.0f * a);
          g2t = lw_dot(rtry, ept) / lw_nrm(ept);
          gnt = sqrt(g1t * g1t + g2t * g2t);
          if (gnt < gn) break;
          lam = 0.5f * lam;
        }
        if (gnt >= gn) {
          conv = gn <= 16.0f * tolf;
          if (!conv && (g1 > 0.0f || fabs(u) >= 1.25f))
            status = LPT_NEWTON_OUTSIDE;
          break;
        }
        al = altry; u = utry; ep = ept; epp = eppt; rr = rtry;
        g1 = g1t; g2 = g2t; gn = gnt;
      }
      if (conv) status = LPT_NEWTON_CONVERGED;
      if (status != LPT_NEWTON_NOT_CONVERGED) break;
    }
    if (status == LPT_NEWTON_NOT_CONVERGED) cnt[LPT_ST_UNRESOLVED] += 1;
    if (!conv) continue;
    if (fabs(u) > 1.0f || al < 0.0f || al > 1.0f || al >= *alpha) continue;
    if (lw_dot(dx, rr) >= 0.0f) continue;
    real s, t;
    lw_edge_st(k, u, &s, &t);
    real3 x, xs, xt, nfe, ns, nt;
    if (!lw_frame(g, f, s, t, &x, &xs, &xt, &nfe, 0, &ns, &nt) ||
        lw_dot(rr, nfe) > tolf)
      continue;
    *alpha = al;
    *found = 1;
    *nc = lw_scale(-1.0f / lw_nrm(rr), rr);
  }
}

/* accept_root */
static inline int lw_accept(real s, real t, real al, real3 n, real3 dx,
                         real tol_st, real tol_al, thread real *alpha,
                         thread int *found, thread real3 *nc) {
  const int ok = fabs(s) <= 1.0f + tol_st && fabs(t) <= 1.0f + tol_st &&
    al >= -tol_al && al <= 1.0f && lw_dot(dx, n) > 0.0f;
  if (ok && lw_max(0.0f, al) < *alpha) {
    *alpha = lw_max(0.0f, al);
    *found = 1;
    *nc = n;
  }
  return ok;
}

/* lpt_wall_facet_intersect, x0g in global coordinates */
static inline void lw_facet_intersect(lw_geom g, int f, real3 x0g, real3 dx,
                                   real a, real sg, real tg, int have_guess,
                                   thread real *alpha, thread int *found,
                                   thread real3 *nc, thread int *is_rim,
                                   real tie, thread uint *cnt) {
  const real eps = LPT_EPS;
  const real half_pi = 1.5707963267948966f;
  *found = 0;
  *is_rim = 0;
  *nc = lw_v3(0.0f, 0.0f, 0.0f);
  const real dxl = lw_nrm(dx);
  if (!(dxl > 0.0f)) return;
  const real3 dxn = lw_scale(1.0f / dxl, dx);
  const real tol_st = 1.0e3f * eps;
  const real tol_al = 64.0f * eps * lw_fr(g, f, LPT_FR_SCALE) / dxl;
  const real3 x0 = lw_sub(x0g, lw_fr3(g, f, LPT_FR_ORIGIN));
  const real3 cl = lw_fr3(g, f, LPT_FR_CENTER);
  const real3 ax0 = lw_fr3(g, f, LPT_FR_AXIS);

  real a_in, a_out;
  if (!lw_seg_sphere(x0, dx, cl, lw_fr(g, f, LPT_FR_RADIUS) + a, &a_in,
                     &a_out) || a_in > *alpha)
    return;
  real hmax = lw_max(lw_dot(lw_sub(x0, cl), ax0),
                     lw_dot(lw_sub(lw_add(x0, dx), cl), ax0));
  if (hmax < lw_fr(g, f, LPT_FR_SLAB) - a) return;

  real al_r = *alpha;
  int fd_r = 0;
  real3 nc_r = lw_v3(0.0f, 0.0f, 0.0f);
  if (a > 0.0f) lw_rim(g, f, x0, dx, a, &al_r, &fd_r, &nc_r, cnt);

  real boxes[4 * LPT_PATCH_STACK];
  int depth[LPT_PATCH_STACK];
  int top = 0;
  boxes[0] = -1.0f; boxes[1] = 1.0f; boxes[2] = -1.0f; boxes[3] = 1.0f;
  depth[0] = 0;
  while (top >= 0) {
    real box[4];
    for (int i = 0; i < 4; i++) box[i] = boxes[4 * top + i];
    const int dep = depth[top];
    top--;

    real3 c, axis;
    real r, theta, slab[2];
    if (dep == 0) {
      c = cl;
      r = lw_fr(g, f, LPT_FR_RADIUS);
      axis = ax0;
      theta = lw_fr(g, f, LPT_FR_ANGLE);
      slab[0] = lw_fr(g, f, LPT_FR_SLAB);
      slab[1] = lw_fr(g, f, LPT_FR_SLAB + 1);
    } else {
      lw_patch_bounds(g, f, box, &c, &r, &axis, &theta, slab);
    }
    if (!lw_seg_sphere(x0, dx, c, r + a, &a_in, &a_out) || a_in > *alpha)
      continue;
    hmax = lw_max(lw_dot(lw_sub(x0, c), axis),
                  lw_dot(lw_sub(lw_add(x0, dx), c), axis));
    if (hmax < slab[0] - a) continue;

    const real cosang = lw_dot(dxn, axis);
    int graph = 0;
    if (theta < half_pi) {
      const real sth = sin(theta);
      if (cosang < -sth) continue;
      graph = cosang > sth;
    }

    int subdivide = 0;
    const real sm = 0.5f * (box[0] + box[1]);
    const real tm = 0.5f * (box[2] + box[3]);
    if (graph || dep == LPT_PATCH_MAX_DEPTH) {
      real s0 = sm, t0 = tm;
      if (dep == 0 && have_guess) {
        s0 = lw_clamp(sg, -1.0f, 1.0f);
        t0 = lw_clamp(tg, -1.0f, 1.0f);
      }
      real3 xg, xs, xt, ng, ns, nt;
      real al0 = 0.5f * (a_in + a_out);
      if (lw_frame(g, f, s0, t0, &xg, &xs, &xt, &ng, 0, &ns, &nt)) {
        const real den = lw_dot(dx, ng);
        if (fabs(den) > eps * dxl)
          al0 = lw_dot(lw_sub(lw_sub(xg, lw_scale(a, ng)), x0), ng) / den;
      }
      real s, t, al;
      real3 nrm;
      int iters;
      int status = lw_newton(g, f, x0, dx, a, box, s0, t0, al0, 0, &s, &t,
                             &al, &nrm, &iters);
      if ((status == LPT_NEWTON_SINGULAR ||
           status == LPT_NEWTON_NOT_CONVERGED) && a > 0.0f &&
          lw_fr(g, f, LPT_FR_PLANAR) == 0.0f) {
        status = lw_newton(g, f, x0, dx, a, box, s0, t0, al0, 1, &s, &t,
                           &al, &nrm, &iters);
      }
      int valid = 0;
      if (status == LPT_NEWTON_CONVERGED)
        valid = lw_accept(s, t, al, nrm, dx, tol_st, tol_al, alpha, found,
                          nc);
      if (!valid && !graph && dep == LPT_PATCH_MAX_DEPTH) {
        status = lw_newton(g, f, x0, dx, a, box, sm, tm, a_out, a > 0.0f,
                           &s, &t, &al, &nrm, &iters);
        if (status == LPT_NEWTON_CONVERGED)
          valid = lw_accept(s, t, al, nrm, dx, tol_st, tol_al, alpha,
                            found, nc);
        if (!valid) cnt[LPT_ST_UNRESOLVED] += 1;
      } else if (!valid) {
        if (status == LPT_NEWTON_CONVERGED ||
            (status == LPT_NEWTON_OUTSIDE && graph)) {
          /* the unique crossing is outside this patch or this step */
        } else if (dep < LPT_PATCH_MAX_DEPTH) {
          subdivide = 1;
        } else {
          cnt[LPT_ST_UNRESOLVED] += 1;
        }
      }
    } else {
      subdivide = 1;
    }

    if (subdivide) {
      for (int k = 0; k < 4; k++) {
        top++;
        depth[top] = dep + 1;
        boxes[4 * top] = (k == 0 || k == 2) ? box[0] : sm;
        boxes[4 * top + 1] = (k == 0 || k == 2) ? sm : box[1];
        boxes[4 * top + 2] = (k < 2) ? box[2] : tm;
        boxes[4 * top + 3] = (k < 2) ? tm : box[3];
      }
    }
  }

  if (fd_r && (!*found || al_r < *alpha - tie)) {
    *alpha = al_r;
    *nc = nc_r;
    *found = 1;
    *is_rim = 1;
  }
}

/* facet_st_from_rst */
static inline void lw_st_from_rst(int facet, real r, real s, real t,
                               thread real *sf, thread real *tf) {
  if (facet == 1 || facet == 2) { *sf = s; *tf = t; }
  else if (facet == 3 || facet == 4) { *sf = r; *tf = t; }
  else { *sf = r; *tf = s; }
  *sf = lw_clamp(*sf, -1.0f, 1.0f);
  *tf = lw_clamp(*tf, -1.0f, 1.0f);
}

/* lpt_wall_collide: el is 1-based; returns the number of hits */
static inline int lw_collide(lw_geom g, int el, real3 p0, thread real3 *p1,
                          real a, real e_n, thread real3 *normals,
                          int has_rst, real r0, real r1, real r2,
                          thread uint *cnt) {
  const real eps = LPT_EPS;
  int n_hits = 0;
  if (el < 1 || el > g.nelv) return 0;
  if (g.el_ptr[el - 1] >= g.el_ptr[el]) return 0;
  if (lw_nrm(lw_sub(*p1, p0)) + a > g.el_radius[el - 1])
    cnt[LPT_ST_REACH] += 1;

  real3 q0 = p0, q1 = *p1;
  for (int bounce = 1; bounce <= LPT_WALL_MAX_BOUNCES + 1; bounce++) {
    const real3 dx = lw_sub(q1, q0);
    if (!(lw_nrm(dx) > 0.0f)) break;
    real alpha = LPT_HUGE;
    int fbest = 0;
    int rim_best = 0;
    real3 nbest = lw_v3(0.0f, 0.0f, 0.0f);
    const real qmax = lw_max(fabs(q0.x), lw_max(fabs(q0.y), fabs(q0.z)));
    const real a_tie = (128.0f * eps * (qmax + lw_nrm(dx)) + 1.0e-3f * a) /
      lw_nrm(dx);
    for (int k = g.el_ptr[el - 1]; k < g.el_ptr[el]; k++) {
      const int f = g.el_facets[k - 1];
      int have_guess = 0;
      real sg = 0.0f, tg = 0.0f;
      if (bounce == 1 && has_rst && g.fi[2 * (f - 1)] == el &&
          g.fi[2 * (f - 1) + 1] > 0) {
        lw_st_from_rst(g.fi[2 * (f - 1) + 1], r0, r1, r2, &sg, &tg);
        have_guess = 1;
      }
      real al_c = rim_best ? alpha + a_tie : alpha;
      int found, rim;
      real3 nc;
      lw_facet_intersect(g, f, q0, dx, a, sg, tg, have_guess, &al_c, &found,
                         &nc, &rim, a_tie, cnt);
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

    const real tol_off = 16.0f * eps * lw_fr(g, fbest, LPT_FR_SCALE);
    if (n_hits == LPT_WALL_MAX_BOUNCES) {
      q1 = lw_sub(lw_add(q0, lw_scale(alpha, dx)), lw_scale(tol_off, nbest));
      break;
    }
    real3 rem = lw_scale(1.0f - alpha, dx);
    q0 = lw_sub(lw_add(q0, lw_scale(alpha, dx)), lw_scale(tol_off, nbest));
    const real rn = lw_dot(rem, nbest);
    if (rn > 0.0f) rem = lw_sub(rem, lw_scale((1.0f + e_n) * rn, nbest));
    q1 = lw_add(q0, rem);
    normals[n_hits] = nbest;
    n_hits++;
  }
  *p1 = q1;
  return n_hits;
}

/* v <- v - (1 + e_n)(v . n) n, on three global arrays at index p */
static inline void lw_reflect(device real *vx, device real *vy, device real *vz,
                           int p, real nx, real ny, real nz, real e_n) {
  const real vn = (1.0f + e_n) * (vx[p] * nx + vy[p] * ny + vz[p] * nz);
  vx[p] -= vn * nx;
  vy[p] -= vn * ny;
  vz[p] -= vn * nz;
}

/** Collision search; writes positions, hit counts and contact normals */
kernel void lpt_wall_collide_kernel(device const float *xf [[ buffer(0) ]],
                                    device const float *fr [[ buffer(1) ]],
                                    device const int *fi [[ buffer(2) ]],
                                    device const float *zg [[ buffer(3) ]],
                                    device const int *el_ptr [[ buffer(4) ]],
                                    device const int *el_facets
                                    [[ buffer(5) ]],
                                    device const float *el_radius
                                    [[ buffer(6) ]],
                                    device const int *el_list [[ buffer(7) ]],
                                    device const float *rst [[ buffer(8) ]],
                                    device const float *x_old [[ buffer(9) ]],
                                    device const float *y_old
                                    [[ buffer(10) ]],
                                    device const float *z_old
                                    [[ buffer(11) ]],
                                    device float *x [[ buffer(12) ]],
                                    device float *y [[ buffer(13) ]],
                                    device float *z [[ buffer(14) ]],
                                    device const float *d [[ buffer(15) ]],
                                    device int *hits [[ buffer(16) ]],
                                    device float *normals [[ buffer(17) ]],
                                    device atomic_uint *stats
                                    [[ buffer(18) ]],
                                    constant int &lx [[ buffer(19) ]],
                                    constant int &nelv [[ buffer(20) ]],
                                    constant int &n [[ buffer(21) ]],
                                    constant float &e_n [[ buffer(22) ]],
                                    constant int &has_rst [[ buffer(23) ]],
                                    uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;
  const int p = (int) idx;
  lw_geom g;
  g.xf = xf; g.fr = fr; g.fi = fi; g.zg = zg; g.el_ptr = el_ptr;
  g.el_facets = el_facets; g.el_radius = el_radius; g.lx = lx; g.nelv = nelv;

  uint cnt[LPT_NSTATS];
  for (int i = 0; i < LPT_NSTATS; i++) cnt[i] = 0;
  hits[p] = 0;
  if (el_list[p] >= 0) {
    real3 p1 = lw_v3(x[p], y[p], z[p]);
    real3 nrm[LPT_WALL_MAX_BOUNCES];
    const int nh = lw_collide(g, el_list[p] + 1,
                              lw_v3(x_old[p], y_old[p], z_old[p]), &p1,
                              0.5f * d[p], e_n, nrm, has_rst,
                              has_rst ? rst[3 * p] : 0.0f,
                              has_rst ? rst[3 * p + 1] : 0.0f,
                              has_rst ? rst[3 * p + 2] : 0.0f, cnt);
    x[p] = p1.x;
    y[p] = p1.y;
    z[p] = p1.z;
    hits[p] = nh;
    for (int k = 0; k < nh; k++) {
      normals[3 * (LPT_WALL_MAX_BOUNCES * p + k)] = nrm[k].x;
      normals[3 * (LPT_WALL_MAX_BOUNCES * p + k) + 1] = nrm[k].y;
      normals[3 * (LPT_WALL_MAX_BOUNCES * p + k) + 2] = nrm[k].z;
    }
  }
  for (int i = 0; i < LPT_NSTATS; i++) {
    if (cnt[i] == 0) continue;
    atomic_fetch_add_explicit(&stats[i], cnt[i], memory_order_relaxed);
  }
}

/** Reflect the velocity, the step's velocity and acceleration, and their
    histories about the stored contact normals */
kernel void lpt_wall_reflect_kernel(device const int *hits [[ buffer(0) ]],
                                    device const float *normals
                                    [[ buffer(1) ]],
                                    device float *u [[ buffer(2) ]],
                                    device float *v [[ buffer(3) ]],
                                    device float *w [[ buffer(4) ]],
                                    device float *u_old [[ buffer(5) ]],
                                    device float *v_old [[ buffer(6) ]],
                                    device float *w_old [[ buffer(7) ]],
                                    device float *acc_x [[ buffer(8) ]],
                                    device float *acc_y [[ buffer(9) ]],
                                    device float *acc_z [[ buffer(10) ]],
                                    device float *u_lag [[ buffer(11) ]],
                                    device float *v_lag [[ buffer(12) ]],
                                    device float *w_lag [[ buffer(13) ]],
                                    device float *acc_xlag [[ buffer(14) ]],
                                    device float *acc_ylag [[ buffer(15) ]],
                                    device float *acc_zlag [[ buffer(16) ]],
                                    device float *u_laglag [[ buffer(17) ]],
                                    device float *v_laglag [[ buffer(18) ]],
                                    device float *w_laglag [[ buffer(19) ]],
                                    device float *acc_xlaglag
                                    [[ buffer(20) ]],
                                    device float *acc_ylaglag
                                    [[ buffer(21) ]],
                                    device float *acc_zlaglag
                                    [[ buffer(22) ]],
                                    constant int &n [[ buffer(23) ]],
                                    constant int &lag_len [[ buffer(24) ]],
                                    constant float &e_n [[ buffer(25) ]],
                                    uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;
  const int p = (int) idx;
  for (int k = 0; k < hits[p]; k++) {
    const int i = 3 * (LPT_WALL_MAX_BOUNCES * p + k);
    const float nx = normals[i], ny = normals[i + 1], nz = normals[i + 2];
    lw_reflect(u, v, w, p, nx, ny, nz, e_n);
    lw_reflect(u_old, v_old, w_old, p, nx, ny, nz, e_n);
    lw_reflect(acc_x, acc_y, acc_z, p, nx, ny, nz, e_n);
    if (lag_len >= 1) {
      lw_reflect(u_lag, v_lag, w_lag, p, nx, ny, nz, e_n);
      lw_reflect(acc_xlag, acc_ylag, acc_zlag, p, nx, ny, nz, e_n);
    }
    if (lag_len >= 2) {
      lw_reflect(u_laglag, v_laglag, w_laglag, p, nx, ny, nz, e_n);
      lw_reflect(acc_xlaglag, acc_ylaglag, acc_zlaglag, p, nx, ny, nz, e_n);
    }
  }
}
