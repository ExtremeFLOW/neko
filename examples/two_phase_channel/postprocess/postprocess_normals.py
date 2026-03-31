"""
postprocess_normals.py — Normal field quality diagnostics for σ=0 CDI runs.

Motivation
----------
The CDI (Compact Diffuse Interface) compression term is ∇·(γ φ(1-φ) n̂).
The unit normal n̂ = ∇φ/|∇φ| is computed from the SEM derivative operators.
If n̂ points in the wrong direction, the sharpening acts away from the interface
even if γ is large — the drop interface will then broaden over time.

In the SEM, the Lagrange endpoint derivative amplification D[N,N] ≈ 14 (N=7)
amplifies any discontinuity between the face-averaged n̂ and the first-interior n̂.
This creates a large intra-element kink in n̂ near element faces. The two figures
below directly measure this effect and its consequence on interface sharpness.

Two figures per run:

  1. normals_angular_dev_<run>.png
     Angular deviation δ = arccos(n̂_computed · n̂_ideal) in degrees at every
     interface GLL point in the x-y midplane, for each of the three snapshots
     (first / mid / last) shown as separate panels.

     n̂_computed = ∇φ/|∇φ| from element-local numpy.gradient (same as SEM but
     on a regular GLL grid — an approximation, but captures the kink pattern).
     n̂_ideal = inward radial direction from the auto-detected drop centroid,
     i.e. -(x−x_c, y−y_c) / r. For a spherical drop, this is exactly correct.

     Healthy run: δ ≈ 0 everywhere.
     SEM kink hypothesis: large δ concentrated at element face crossings —
     visible as a regular striped pattern aligned with the element grid.
     Mean δ ~ 60° → normals are severely wrong; CDI sharpening will fail.

  2. normals_phi_profile_<run>.png
     1D φ profile through the top of the drop (vertical cut at x = x_centroid)
     for first / mid / last snapshot, all shifted to s = y − y_interface.
     Ideal tanh profile plotted for reference. 10-90% width measured and printed.
     Expected width = 4ε · arctanh(0.8) ≈ 4.4ε.

     Profile stays unchanged → CDI sharpening is working (γ adequate, n̂ adequate).
     Profile broadens over time → CDI sharpening failing, either because γ is too
     low OR because n̂ is wrong (both causes produce the same observable symptom).
     Combined with Figure 1: if δ is large AND the profile broadens, bad n̂ is
     the primary cause. If δ is small but the profile still broadens, γ is too low.

Usage:
    cd examples/two_phase_channel
    source ../../setup-env-channel.sh --egidius
    python3 postprocess/postprocess_normals.py --run channel_p3_sigma0 --R 0.4 --eps 0.04 --mesh l2
    python3 postprocess/postprocess_normals.py --run channel_l3_sigma0 --R 0.4 --eps 0.03 --mesh l3
    python3 postprocess/postprocess_normals.py --run channel_p2_sigma0_eps053 --R 0.4 --eps 0.053 --mesh p2
"""

import argparse
import glob
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.msh import Mesh as msh_c
from pysemtools.datatypes.field import Field as field_c

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--run',  default='channel_p3_sigma0')
parser.add_argument('--R',   type=float, default=0.4)
parser.add_argument('--eps', type=float, default=0.04)
parser.add_argument('--mesh', choices=['p1', 'p2', 'l2', 'l3', 'l4'], default='l2')
parser.add_argument('--grad-thresh', type=float, default=None,
                    help='|∇φ| threshold to identify interface GLL points '
                         '(default: 1/(8ε), half the peak gradient of ideal tanh)')
parser.add_argument('--dpi', type=int, default=160)
args = parser.parse_args()

SIM_DIR = '/lscratch/sieburgh/simulations'
RUN_DIR = f'{SIM_DIR}/{args.run}'
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)
comm = MPI.COMM_WORLD

R   = args.R
EPS = args.eps
LLX = 4.0 * math.pi
LLZ = 4.0 / 3.0 * math.pi
Z_C = LLZ / 2.0

_NZ = {'p1': 27, 'p2': 36, 'l2': 48, 'l3': 64, 'l4': 96}
NZ_ELEMS = _NZ[args.mesh]
DXZ = LLZ / NZ_ELEMS   # approximate element size

_SPINUP_RUNS = {
    'p1': 'channel_single_phase',
    'p2': 'channel_p2_single_phase',
    'l2': 'channel_p3_single_phase',
    'l3': 'channel_l3_single_phase',
    'l4': 'channel_l4_single_phase',
}

# Gradient threshold: 1/(8ε) = half the peak |∇φ| of the ideal tanh profile.
# Points above this threshold are in the interface band.
GRAD_THRESH = args.grad_thresh if args.grad_thresh else 1.0 / (8.0 * EPS)

# ---------------------------------------------------------------------------
# Load mesh
# ---------------------------------------------------------------------------
_spinup_f0 = os.path.join(SIM_DIR, _SPINUP_RUNS[args.mesh], 'field0.f00000')
_candidates = [os.path.join(RUN_DIR, 'field0.f00000'), _spinup_f0]
mesh_file = next((f for f in _candidates if os.path.exists(f)), None)
if mesh_file is None:
    raise FileNotFoundError('No mesh file found. Check --mesh flag and spin-up directory.')
print(f'Reading mesh from: {mesh_file}')

xyz_data = preadnek(mesh_file, comm)
msh = msh_c(comm, data=xyz_data)
del xyz_data

lz_ = msh.x.shape[1]
ly_ = msh.x.shape[2]
lx_ = msh.x.shape[3]
iz_mid = lz_ // 2
iy_mid = ly_ // 2

# X-Y midplane slice: select elements whose z-centroid is near LLZ/2
z_ec    = msh.z[:, iz_mid, iy_mid, lx_ // 2]
mask_xy = np.abs(z_ec - Z_C) < 0.6 * DXZ
n_el_sl = mask_xy.sum()
print(f'x-y slice: {n_el_sl} elements')

# Flattened GLL coordinates of the x-y slice
x_gll = msh.x[mask_xy, iz_mid, :, :].reshape(-1)   # (n_el_sl * ly_ * lx_,)
y_gll = msh.y[mask_xy, iz_mid, :, :].reshape(-1)

# Element centroid coords (for spatial lookups)
x_ec_sl = msh.x[mask_xy, iz_mid, iy_mid, lx_ // 2]
y_ec_sl = msh.y[mask_xy, iz_mid, iy_mid, lx_ // 2]

# Element boundary lines (for figure 1)
x_lines = np.unique(np.round(msh.x[mask_xy, iz_mid, :, 0].reshape(-1), 5))
y_lines = np.unique(np.round(msh.y[mask_xy, iz_mid, 0, :].reshape(-1), 5))

# Triangulation (for φ=0.5 contour in figure 1)
def _build_triang(xs, ys):
    n_el, nj, ni = xs.shape
    pts_x, pts_y, tris = [], [], []
    offset = 0
    for e in range(n_el):
        pts_x.append(xs[e].flatten())
        pts_y.append(ys[e].flatten())
        for j in range(nj - 1):
            for i in range(ni - 1):
                p0 = offset + j * ni + i
                tris.extend([[p0, p0+1, offset+(j+1)*ni+i+1],
                              [p0, offset+(j+1)*ni+i+1, offset+(j+1)*ni+i]])
        offset += nj * ni
    return tri.Triangulation(np.concatenate(pts_x), np.concatenate(pts_y),
                             np.array(tris))

triang_xy = _build_triang(msh.x[mask_xy, iz_mid, :, :],
                           msh.y[mask_xy, iz_mid, :, :])

# Copy slice element coordinates for gradient computation, then free full mesh
xs_sl = msh.x[mask_xy].copy()   # (n_el_sl, lz_, ly_, lx_)
ys_sl = msh.y[mask_xy].copy()
zs_sl = msh.z[mask_xy].copy()
del msh

# ---------------------------------------------------------------------------
# Gradient helper: element-local ∇φ (gx, gy) in the x-y midplane
# ---------------------------------------------------------------------------
def gradients_2d(phi_arr):
    """Element-local 2-D gradients (gx, gy) at iz_mid for the x-y slice.

    Returns flattened arrays of shape (n_el_sl * ly_ * lx_,) aligned with
    x_gll / y_gll.
    """
    sub = phi_arr[mask_xy]   # (n_el_sl, lz_, ly_, lx_)
    gx  = np.zeros(n_el_sl * ly_ * lx_)
    gy  = np.zeros(n_el_sl * ly_ * lx_)
    for e in range(n_el_sl):
        x1d = xs_sl[e, 0, 0, :]
        y1d = ys_sl[e, 0, :, 0]
        z1d = zs_sl[e, :, 0, 0]
        gz_e, gy_e, gx_e = np.gradient(sub[e], z1d, y1d, x1d)
        s = e * ly_ * lx_
        gx[s:s + ly_ * lx_] = gx_e[iz_mid].reshape(-1)
        gy[s:s + ly_ * lx_] = gy_e[iz_mid].reshape(-1)
    return gx, gy

# ---------------------------------------------------------------------------
# Snapshot selection: first / mid / last
# ---------------------------------------------------------------------------
field_files = sorted(glob.glob(os.path.join(RUN_DIR, 'field0.f[0-9]*')))
if not field_files:
    raise FileNotFoundError(f'No field0.f* files in {RUN_DIR}')

n_files = len(field_files)
indices = sorted(set([0, n_files // 2, n_files - 1]))
snap_files = [field_files[i] for i in indices]
print(f'Loading {len(snap_files)} snapshots ({n_files} total, using first/mid/last) ...')

snaps = []
for fpath in snap_files:
    print(f'  {os.path.basename(fpath)} ...', end=' ', flush=True)
    data    = preadnek(fpath, comm)
    fld     = field_c(comm, data=data)
    t       = fld.t
    phi_arr = fld.fields['scal'][0]
    phi_xy  = phi_arr[mask_xy, iz_mid, :, :].reshape(-1)
    gx, gy  = gradients_2d(phi_arr)
    del data, fld
    snaps.append({'t': t, 'phi_xy': phi_xy, 'gx': gx, 'gy': gy})
    print(f't={t:.3f}')

# ---------------------------------------------------------------------------
# Figure 1: Angular deviation map — one panel per snapshot (first / mid / last)
#
# Why: directly measures the geometric error in n̂ at each interface GLL point.
# n̂_computed = ∇φ/|∇φ| (element-local gradient via numpy.gradient).
# n̂_ideal    = inward radial from drop centroid = -(x−xc, y−yc)/r.
# δ = arccos(n̂_computed · n̂_ideal): 0° = perfect, 90° = perpendicular, 180° = inverted.
# A mean δ ~ 60° means normals are severely wrong → CDI sharpening term
# ∇·(γ φ(1-φ) n̂) is not pointing toward the interface → interface broadens.
# SEM kink signature: large δ in stripes aligned with element face crossings.
# ---------------------------------------------------------------------------
print('\nBuilding Figure 1: angular deviation map (all 3 snapshots) ...')


def _compute_deviation(snap):
    """Return (x_valid, y_valid, delta, stats_dict) for one snapshot."""
    phi_xy = snap['phi_xy']
    gx     = snap['gx']
    gy     = snap['gy']

    # Auto-detect drop centroid from φ > 0.5 weighted average
    _in_drop = phi_xy > 0.5
    if _in_drop.sum() > 10:
        _w  = phi_xy[_in_drop] - 0.5
        x_c = float(np.average(x_gll[_in_drop], weights=_w))
        y_c = float(np.average(y_gll[_in_drop], weights=_w))
    else:
        x_c, y_c = LLX / 2.0, 0.0

    grad_mag = np.sqrt(gx**2 + gy**2)
    in_iface = (grad_mag > GRAD_THRESH) & (phi_xy > 0.05) & (phi_xy < 0.95)
    dx = x_gll - x_c
    dy = y_gll - y_c
    r  = np.sqrt(dx**2 + dy**2)
    valid = in_iface & (r > R * 0.25)

    nx_c = gx[valid] / grad_mag[valid]
    ny_c = gy[valid] / grad_mag[valid]
    nx_i = -dx[valid] / r[valid]
    ny_i = -dy[valid] / r[valid]

    dot   = np.clip(nx_c * nx_i + ny_c * ny_i, -1.0, 1.0)
    delta = np.degrees(np.arccos(dot))
    stats = dict(mean=delta.mean(), median=float(np.median(delta)),
                 p90=float(np.percentile(delta, 90)), mx=delta.max(),
                 n=valid.sum(), x_c=x_c, y_c=y_c)
    return x_gll[valid], y_gll[valid], phi_xy, delta, stats


dev_results = []
for snap in snaps:
    xv, yv, phi_xy, delta, stats = _compute_deviation(snap)
    dev_results.append((snap['t'], xv, yv, phi_xy, delta, stats))
    print(f'  t={snap["t"]:.2f}: centroid=({stats["x_c"]:.3f},{stats["y_c"]:.3f})  '
          f'n={stats["n"]}  δ mean={stats["mean"]:.1f}°  median={stats["median"]:.1f}°  '
          f'90th={stats["p90"]:.1f}°  max={stats["mx"]:.1f}°')

# Shared colour scale: 98th-percentile of ALL snapshots combined
all_deltas = np.concatenate([d for _, _, _, _, d, _ in dev_results])
vmax_sc = float(np.percentile(all_deltas, 98))

# All three snapshots on the same panel — the drop advects ~1 unit/TU so the
# three clusters are well-separated in x and can coexist without overlap.
fig, ax = plt.subplots(figsize=(13, 2.8), dpi=args.dpi)
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.04, right=0.87, top=0.80, bottom=0.10)

contour_colors = ['#1565C0', '#EF6C00', '#B71C1C']   # blue / orange / red: first/mid/last

for (t, xv, yv, phi_xy, delta, stats), cc in zip(dev_results, contour_colors):
    ax.scatter(xv, yv, c=delta, cmap='hot_r', s=3, vmin=0.0, vmax=vmax_sc,
               rasterized=True, zorder=3, alpha=0.9)
    ax.tricontour(triang_xy, phi_xy, levels=[0.5],
                  colors=[cc], linewidths=1.5, zorder=4)
    # Annotate centroid with time label
    ax.text(stats['x_c'], stats['y_c'] + R + 0.08,
            f't={t:.1f}', fontsize=7, ha='center', color=cc, fontweight='bold')

for xline in x_lines:
    ax.axvline(xline, color='#666666', lw=0.12, alpha=0.25, zorder=1)
for yline in y_lines:
    ax.axhline(yline, color='#666666', lw=0.12, alpha=0.25, zorder=1)

ax.set_xlim(0, LLX)
ax.set_ylim(-1.0, 1.0)
ax.set_xlabel('x', fontsize=8)
ax.set_ylabel('y', fontsize=8)
ax.tick_params(labelsize=7)

sm = ScalarMappable(cmap='hot_r', norm=Normalize(vmin=0.0, vmax=vmax_sc))
sm.set_array([])
cbar_ax = fig.add_axes([0.88, 0.10, 0.015, 0.68])
cb = fig.colorbar(sm, cax=cbar_ax)
cb.set_label('δ  [°]', fontsize=9)
cb.ax.tick_params(labelsize=8)

# Per-snapshot stats in title
stats_lines = '  |  '.join(
    f't={t:.1f}: mean={s["mean"]:.0f}°  p90={s["p90"]:.0f}°  max={s["mx"]:.0f}°'
    for t, _, _, _, _, s in dev_results
)
ax.set_title(
    f'Angular deviation  δ = arccos(n̂_computed · n̂_ideal)   —   {args.run}\n'
    f'n̂_computed = ∇φ/|∇φ|.  n̂_ideal = inward radial from centroid.  '
    f'Colored contours: φ=0.5 at each time.\n'
    + stats_lines,
    fontsize=7.5, fontweight='bold', loc='left')

out1 = os.path.join(OUT_DIR, f'normals_angular_dev_{args.run}.png')
plt.savefig(out1, dpi=args.dpi, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out1}')

# ---------------------------------------------------------------------------
# Figure 2: φ interface profile width over time
# ---------------------------------------------------------------------------
print('\nBuilding Figure 2: φ profile width ...')

# Expected 10-90% width for ideal tanh: w = 4ε × arctanh(0.8)
width_ideal = 4.0 * EPS * float(np.arctanh(0.8))
print(f'  Expected 10-90% width: {width_ideal:.4f}  '
      f'= {width_ideal/EPS:.3f} ε  (ideal tanh)')

# Ideal profile for plotting (top of drop: φ decreases going up, s = y - y_int)
s_ref   = np.linspace(-6.0 * EPS, 6.0 * EPS, 600)
phi_ref = 0.5 * (1.0 - np.tanh(s_ref / (2.0 * EPS)))

fig, ax = plt.subplots(figsize=(7, 5), dpi=args.dpi)
fig.patch.set_facecolor('white')
colors = ['#1565C0', '#EF6C00', '#B71C1C']   # dark blue / orange / dark red

for snap, color in zip(snaps, colors):
    phi_s = snap['phi_xy']

    # Centroid for this snapshot
    _in = phi_s > 0.5
    if _in.sum() > 10:
        _w  = phi_s[_in] - 0.5
        xc  = float(np.average(x_gll[_in], weights=_w))
        yc  = float(np.average(y_gll[_in], weights=_w))
    else:
        xc, yc = LLX / 2.0, 0.0

    # Select GLL points near x = xc (within ±DXZ) in the upper half of the drop
    near_x   = np.abs(x_gll - xc) < DXZ
    upper    = y_gll > yc + R * 0.2
    mask_cut = near_x & upper

    if mask_cut.sum() < 8:
        print(f'  t={snap["t"]:.2f}: too few points near x_c={xc:.2f}, skipping')
        continue

    y_cut   = y_gll[mask_cut]
    phi_cut = phi_s[mask_cut]
    order   = np.argsort(y_cut)
    y_cut   = y_cut[order]
    phi_cut = phi_cut[order]

    # Find interface position (φ = 0.5, going from φ>0.5 below to φ<0.5 above)
    sign_cross = np.where(np.diff(np.sign(phi_cut - 0.5)) < 0)[0]
    if len(sign_cross) == 0:
        print(f'  t={snap["t"]:.2f}: no φ=0.5 crossing found')
        continue
    ci    = sign_cross[-1]   # topmost crossing (upper interface)
    y0, y1 = y_cut[ci], y_cut[ci + 1]
    p0, p1 = phi_cut[ci], phi_cut[ci + 1]
    y_int = y0 + (0.5 - p0) * (y1 - y0) / (p1 - p0)

    s_cut = y_cut - y_int   # shift: s=0 at interface, s>0 outside

    # Measure 10-90% width
    def _interp_crossing(phi_arr, y_arr, level):
        """Find y where phi crosses level (decreasing)."""
        idx = np.where(np.diff(np.sign(phi_arr - level)) < 0)[0]
        if len(idx) == 0:
            return None
        i = idx[-1]
        return y_arr[i] + (level - phi_arr[i]) * (y_arr[i+1] - y_arr[i]) / (phi_arr[i+1] - phi_arr[i])

    y09 = _interp_crossing(phi_cut, y_cut, 0.9)
    y01 = _interp_crossing(phi_cut, y_cut, 0.1)
    if y09 is not None and y01 is not None:
        width = y01 - y09
        wstr  = f'  w = {width:.4f} = {width/EPS:.2f}ε'
    else:
        wstr = ''

    print(f'  t={snap["t"]:.2f}: y_int={y_int:.3f}{wstr}')
    ax.plot(s_cut / EPS, phi_cut, color=color, lw=1.3, alpha=0.85,
            label=f't = {snap["t"]:.2f} TU{wstr}')

# Reference ideal profile
ax.plot(s_ref / EPS, phi_ref, 'k--', lw=1.6,
        label=f'Ideal tanh   w = {width_ideal/EPS:.2f}ε = {width_ideal:.4f}')

# Guide lines
ax.axvline(0,   color='k',    lw=0.6, ls=':')
ax.axhline(0.5, color='k',    lw=0.6, ls=':')
ax.axhline(0.1, color='gray', lw=0.5, ls=':', alpha=0.6)
ax.axhline(0.9, color='gray', lw=0.5, ls=':', alpha=0.6)
ax.text(5.5, 0.12, '10%', fontsize=7, color='gray', va='bottom')
ax.text(5.5, 0.91, '90%', fontsize=7, color='gray', va='bottom')

ax.set_xlim(-6, 6)
ax.set_ylim(-0.05, 1.05)
ax.set_xlabel('s / ε   (s = y − y_interface,  s > 0 outside drop)', fontsize=10)
ax.set_ylabel('φ', fontsize=10)
ax.legend(fontsize=8, loc='lower left')
ax.grid(True, alpha=0.3)
ax.set_title(
    f'φ interface profile — {args.run}   (R={R}, ε={EPS})\n'
    f'Vertical cut at x ≈ x_centroid through top of drop.  '
    f'Broadening → CDI sharpening failing (low γ and/or bad n̂).',
    fontsize=9, fontweight='bold')

out2 = os.path.join(OUT_DIR, f'normals_phi_profile_{args.run}.png')
plt.savefig(out2, dpi=args.dpi, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out2}')

print('\nAll done.')
