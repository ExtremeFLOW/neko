"""
postprocess_sigma0.py — Diagnostic figures for σ=0 CDI runs.

Produces:
  1. Time-series figure  (reads ekin.csv — always fast)
        κ_rms(t) with 2/R reference line
        φ extrema(t)
        E_kin(t) — turbulence stability check

  2. Field snapshots  (first / mid / last field files)
        Row 0: φ             — where is the drop
        Row 1: |u|           — turbulent flow context
        Row 2: n̂_y          — normal field smoothness (kinks visible at element faces)
        Row 3: κ             — curvature artifact

  3. Zoom + quiver figure  (--normals flag, last snapshot only)
        ~6×4 element patch around top of drop
        Left:  φ background + n̂ quiver (only where |∇φ| > threshold)
        Right: κ background + n̂ quiver
        Element boundaries drawn thick — kink locations clearly marked

Usage:
    python3 postprocess_sigma0.py --run channel_p3_sigma0 --R 0.4 --eps 0.04 --mesh l2
    python3 postprocess_sigma0.py --run channel_l3_sigma0 --R 0.4 --eps 0.03 --mesh l3
    python3 postprocess_sigma0.py --run channel_p3_sigma0 --R 0.4 --eps 0.04 --mesh l2 --normals
    python3 postprocess_sigma0.py --run channel_p2_sigma0_eps053 --R 0.4 --eps 0.053 --mesh p2 --no-snapshots
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
import pandas as pd
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
parser.add_argument('--run', default='channel_test_sigma0_long',
                    help='Run directory name under /lscratch/sieburgh/simulations/')
parser.add_argument('--R', type=float, default=0.4,
                    help='Drop radius (default: 0.4)')
parser.add_argument('--eps', type=float, default=0.07,
                    help='Interface width ε')
parser.add_argument('--drop-center-y', type=float, default=0.0,
                    help='Wall-normal drop centre offset (default: 0.0)')
parser.add_argument('--mesh', choices=['p1', 'p2', 'l2', 'l3', 'l4'], default='p1',
                    help='Mesh preset for z-slice selection '
                         '(p1=81x18x27 nz=27, p2=108x18x36 nz=36, '
                         'l2=144x24x48 nz=48, l3=192x32x64 nz=64, '
                         'l4=288x48x96 nz=96)')
parser.add_argument('--kappa-scale', type=float, default=None,
                    help='Symmetric κ colour limit for snapshots (auto if unset)')
parser.add_argument('--no-snapshots', action='store_true',
                    help='Skip field file loading; produce only the time-series figure')
parser.add_argument('--normals', action='store_true',
                    help='Also produce zoom+quiver figure for normal field analysis')
parser.add_argument('--dpi', type=int, default=160)
args = parser.parse_args()

RUN_DIR = f'/lscratch/sieburgh/simulations/{args.run}'
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
comm    = MPI.COMM_WORLD

R            = args.R
EPS          = args.eps
KAPPA_SPHERE = 2.0 / R
DROP_Y       = args.drop_center_y

LLX   = 4.0 * math.pi
LLZ   = 4.0 / 3.0 * math.pi
X_C   = LLX / 2.0
Z_C   = LLZ / 2.0
_NZ   = {'p1': 27, 'p2': 36, 'l2': 48, 'l3': 64, 'l4': 96}
NZ_ELEMS = _NZ[args.mesh]

run_label = args.run.replace('channel_test_', '').replace('channel_', '').replace('_', ' ')

# ---------------------------------------------------------------------------
# 1. Time-series from ekin.csv
# ---------------------------------------------------------------------------
ekin_path = os.path.join(RUN_DIR, 'ekin.csv')
if not os.path.exists(ekin_path):
    raise FileNotFoundError(f'ekin.csv not found in {RUN_DIR}')

df = pd.read_csv(
    ekin_path, comment='#', header=None,
    names=['t', 'Ekin', 'enst', 'u_max',
           'kappa_max', 'kappa_min', 'kappa_rms',
           'Fst_max', 'phi_min', 'phi_max', 'extra'],
)
print(f'ekin.csv: {len(df)} rows,  t = {df.t.min():.3f} → {df.t.max():.3f}')

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), dpi=args.dpi)
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.10, right=0.97, top=0.93, bottom=0.07, hspace=0.45)

ax1.plot(df.t, df.kappa_rms, color='steelblue', lw=1.5)
ax1.axhline(KAPPA_SPHERE, color='crimson', ls='--', lw=1.2,
            label=f'2/R = {KAPPA_SPHERE:.2f}  (sphere)')
ax1.set_ylabel('κ_rms', fontsize=11)
ax1.set_xlabel('t  [TU]', fontsize=10)
ax1.set_title('Curvature rms', fontsize=11, fontweight='bold')
ax1.legend(fontsize=9)
ax1.tick_params(labelsize=9)
ax1.set_xlim(df.t.min(), df.t.max())

ax2.plot(df.t, df.phi_max, color='steelblue',  lw=1.5, label='φ_max')
ax2.plot(df.t, df.phi_min, color='darkorange', lw=1.5, label='φ_min')
ax2.axhline(1.0, color='gray', ls=':', lw=0.8)
ax2.axhline(0.0, color='gray', ls=':', lw=0.8)
ax2.set_ylabel('φ', fontsize=11)
ax2.set_xlabel('t  [TU]', fontsize=10)
ax2.set_title('Phase field extrema', fontsize=11, fontweight='bold')
ax2.legend(fontsize=9)
ax2.tick_params(labelsize=9)
ax2.set_xlim(df.t.min(), df.t.max())

ax3.plot(df.t, df.Ekin, color='steelblue', lw=1.5)
ax3.set_ylabel('E_kin', fontsize=11)
ax3.set_xlabel('t  [TU]', fontsize=10)
ax3.set_title('Kinetic energy  (turbulence stability)', fontsize=11, fontweight='bold')
ax3.tick_params(labelsize=9)
ax3.set_xlim(df.t.min(), df.t.max())

fig.suptitle(f'{run_label}   R={R}  ε={EPS}', fontsize=12, fontweight='bold')

out_ts = os.path.join(OUT_DIR, f'sigma0_timeseries_{args.run}.png')
plt.savefig(out_ts, dpi=args.dpi, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out_ts}')

# ---------------------------------------------------------------------------
# 2. Field snapshots (optional)
# ---------------------------------------------------------------------------
if args.no_snapshots:
    print('Skipping field snapshots (--no-snapshots).')
    raise SystemExit(0)

field_files = sorted(glob.glob(os.path.join(RUN_DIR, 'field0.f[0-9]*')))
if not field_files:
    print('No field files found — skipping snapshots.')
    raise SystemExit(0)

n = len(field_files)
snap_files = [field_files[i] for i in sorted({0, n // 2, n - 1})]
print(f'\nLoading {len(snap_files)} snapshots for field plots ...')

# --- Mesh ---
print('Reading mesh ...')
xyz_data = preadnek(snap_files[0], comm)
msh = msh_c(comm, data=xyz_data)

lz_ = msh.x.shape[1]
ly_ = msh.x.shape[2]
lx_ = msh.x.shape[3]
iz_mid = lz_ // 2
iy_mid = ly_ // 2

dz_approx = LLZ / NZ_ELEMS
z_ec      = msh.z[:, iz_mid, iy_mid, lx_ // 2]
mask_xy   = np.abs(z_ec - Z_C) < 0.6 * dz_approx
print(f'  x-y slice (z≈{Z_C:.2f}): {mask_xy.sum()} elements')

# --- Triangulation helper ---
def build_triang(xs, ys):
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

triang_xy = build_triang(msh.x[mask_xy, iz_mid, :, :],
                          msh.y[mask_xy, iz_mid, :, :])

x_lines = np.unique(np.round(msh.x[mask_xy, iz_mid, :, 0].reshape(-1), 5))
y_lines = np.unique(np.round(msh.y[mask_xy, iz_mid, 0, :].reshape(-1), 5))

# --- Slice computation helpers ---
def _element_gradients(phi_arr):
    """Element-local ∇φ on x-y slice. Returns (gx, gy) each shape (n_el, ly_, lx_)."""
    sub = phi_arr[mask_xy]
    xs  = msh.x[mask_xy]; ys = msh.y[mask_xy]; zs = msh.z[mask_xy]
    n_el = sub.shape[0]
    gx_out = np.zeros((n_el, ly_, lx_))
    gy_out = np.zeros((n_el, ly_, lx_))
    for e in range(n_el):
        x1d = xs[e, 0, 0, :]; y1d = ys[e, 0, :, 0]; z1d = zs[e, :, 0, 0]
        gz, gy, gx = np.gradient(sub[e], z1d, y1d, x1d)
        gx_out[e] = gx[iz_mid]
        gy_out[e] = gy[iz_mid]
    return gx_out, gy_out

def compute_kappa_slice(phi_arr):
    sub = phi_arr[mask_xy]
    xs  = msh.x[mask_xy]; ys = msh.y[mask_xy]; zs = msh.z[mask_xy]
    n_el = sub.shape[0]
    kap = np.zeros((n_el, ly_, lx_))
    for e in range(n_el):
        x1d = xs[e, 0, 0, :]; y1d = ys[e, 0, :, 0]; z1d = zs[e, :, 0, 0]
        gz, gy, gx = np.gradient(sub[e], z1d, y1d, x1d)
        mag = np.sqrt(gx**2 + gy**2 + gz**2)
        mag = np.where(mag < 1e-3, 1e-3, mag)
        nx_ = gx / mag; ny_ = gy / mag; nz_ = gz / mag
        _, _, dnx = np.gradient(nx_, z1d, y1d, x1d)
        _, dny, _ = np.gradient(ny_, z1d, y1d, x1d)
        dnz, _, _ = np.gradient(nz_, z1d, y1d, x1d)
        kap[e] = -(dnx + dny + dnz)[iz_mid]
    return kap.reshape(-1)

def compute_nhat_y_slice(phi_arr, grad_threshold=0.5):
    """n̂_y on x-y slice. Zero where |∇φ| < threshold (away from interface)."""
    gx, gy = _element_gradients(phi_arr)
    mag = np.sqrt(gx**2 + gy**2)
    nhy = np.where(mag > grad_threshold, gy / np.where(mag > grad_threshold, mag, 1.0), 0.0)
    return nhy.reshape(-1)

def compute_umag_slice(fld):
    """Velocity magnitude |u| on x-y slice. Returns None if vel not in field."""
    try:
        u_arr = fld.fields['vel'][0]
        v_arr = fld.fields['vel'][1]
        w_arr = fld.fields['vel'][2]
    except (KeyError, IndexError):
        return None
    u_xy = u_arr[mask_xy, iz_mid, :, :].reshape(-1)
    v_xy = v_arr[mask_xy, iz_mid, :, :].reshape(-1)
    w_xy = w_arr[mask_xy, iz_mid, :, :].reshape(-1)
    return np.sqrt(u_xy**2 + v_xy**2 + w_xy**2)

# --- Load snapshots ---
snaps = []
for fpath in snap_files:
    print(f'  {os.path.basename(fpath)} ...', end=' ', flush=True)
    data = preadnek(fpath, comm)
    fld  = field_c(comm, data=data)
    t    = fld.t
    phi_arr = fld.fields['scal'][0]
    phi_xy  = phi_arr[mask_xy, iz_mid, :, :].reshape(-1)
    kap_xy  = compute_kappa_slice(phi_arr)
    nhy_xy  = compute_nhat_y_slice(phi_arr)
    umag_xy = compute_umag_slice(fld)
    row = df.iloc[(df.t - t).abs().argsort().iloc[0]]
    print(f't={t:.3f}  κ_rms={row.kappa_rms:.1f}  φ_max={row.phi_max:.4f}')
    snaps.append({'t': t, 'phi_xy': phi_xy, 'kap_xy': kap_xy,
                  'nhy_xy': nhy_xy, 'umag_xy': umag_xy,
                  'phi_arr': phi_arr, 'kappa_rms': row.kappa_rms})

# --- Colour scales ---
kap_absmax = max(np.percentile(np.abs(s['kap_xy']), 99) for s in snaps)
kap_scale  = args.kappa_scale if args.kappa_scale else round(min(kap_absmax, 50.0))

umag_vals = [s['umag_xy'] for s in snaps if s['umag_xy'] is not None]
umag_max  = float(np.percentile(np.concatenate(umag_vals), 99)) if umag_vals else 1.5
umag_max  = round(umag_max, 1)

print(f'\nκ colour limit: ±{kap_scale:.0f}   |u| colour limit: {umag_max:.1f}')

# --- Render snapshots figure: 4 rows × N_snaps cols ---
LEVELS = 100
N_snaps = len(snaps)
n_rows = 4 if umag_vals else 3   # drop |u| row if velocity not available

fig, axes = plt.subplots(n_rows, N_snaps, figsize=(5 * N_snaps, 3.5 * n_rows),
                          dpi=args.dpi)
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.05, right=0.93, top=0.93, bottom=0.05,
                    hspace=0.4, wspace=0.28)
if N_snaps == 1:
    axes = axes.reshape(n_rows, 1)

def _add_panel(ax, triang, vals, cmap, vmin, vmax, title, cb_label, cb_ticks,
               elem_lines=False, phi_vals=None):
    ax.tricontourf(triang, vals, levels=LEVELS, cmap=cmap,
                   vmin=vmin, vmax=vmax, extend='both')
    sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02, aspect=14, fraction=0.03)
    cb.set_label(cb_label, fontsize=8)
    cb.set_ticks(cb_ticks)
    cb.ax.tick_params(labelsize=7)
    if elem_lines:
        for xv in x_lines:
            ax.axvline(xv, color='gray', lw=0.2, alpha=0.35, zorder=1)
        for yv in y_lines:
            ax.axhline(yv, color='gray', lw=0.2, alpha=0.35, zorder=1)
    if phi_vals is not None:
        ax.tricontour(triang, phi_vals, levels=[0.5],
                      colors='limegreen', linewidths=1.2, zorder=2)
    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.set_xlabel('x', fontsize=8)
    ax.set_ylabel('y', fontsize=8)
    ax.tick_params(labelsize=7)

row_phi  = 0
row_umag = 1 if umag_vals else None
row_nhy  = 2 if umag_vals else 1
row_kap  = 3 if umag_vals else 2

for col, s in enumerate(snaps):
    phi_xy  = s['phi_xy']
    kap_xy  = s['kap_xy']
    nhy_xy  = s['nhy_xy']
    umag_xy = s['umag_xy']

    # Row 0: φ
    _add_panel(axes[row_phi, col], triang_xy, phi_xy, 'RdBu_r', 0.0, 1.0,
               f'φ   t = {s["t"]:.2f} TU', 'φ', [0, 0.5, 1.0],
               elem_lines=True)

    # Row 1: |u|  (physical context — no element lines, drop contour for orientation)
    if row_umag is not None and umag_xy is not None:
        _add_panel(axes[row_umag, col], triang_xy, umag_xy, 'plasma',
                   0.0, umag_max,
                   f'|u|   t = {s["t"]:.2f} TU', '|u|',
                   [0, round(umag_max/2, 1), umag_max],
                   elem_lines=False, phi_vals=phi_xy)

    # Row 2: n̂_y  (normal field — element lines key to seeing kinks)
    _add_panel(axes[row_nhy, col], triang_xy, nhy_xy, 'RdBu_r', -1.0, 1.0,
               f'n̂_y   t = {s["t"]:.2f} TU', 'n̂_y', [-1, 0, 1],
               elem_lines=True, phi_vals=phi_xy)

    # Row 3: κ
    _add_panel(axes[row_kap, col], triang_xy, kap_xy, 'seismic',
               -kap_scale, kap_scale,
               f'κ   κ_rms = {s["kappa_rms"]:.1f}  (ref {KAPPA_SPHERE:.1f})',
               'κ', [-kap_scale, 0, kap_scale],
               elem_lines=True, phi_vals=phi_xy)

fig.suptitle(
    f'{run_label}   R={R}  ε={EPS}   snapshots: first / mid / last',
    fontsize=11, fontweight='bold')

out_snaps = os.path.join(OUT_DIR, f'sigma0_snapshots_{args.run}.png')
plt.savefig(out_snaps, dpi=args.dpi, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out_snaps}')

# ---------------------------------------------------------------------------
# 3. Zoom + quiver figure  (--normals flag, last snapshot)
# ---------------------------------------------------------------------------
if not args.normals:
    print('Done.  (use --normals for zoom+quiver figure)')
    raise SystemExit(0)

print('\nBuilding zoom + quiver figure (last snapshot) ...')

last = snaps[-1]
phi_arr = last['phi_arr']

# Element size (approximate, isotropic)
Dxz = LLZ / NZ_ELEMS

# Zoom region: ~6 elements wide × ~4 elements tall, centred on top of drop
x_zoom_c = X_C
y_zoom_c = DROP_Y + R
half_x = 3.0 * Dxz
half_y = 2.0 * Dxz

x_ec = msh.x[:, iz_mid, iy_mid, lx_ // 2]
y_ec = msh.y[:, iz_mid, iy_mid, lx_ // 2]
mask_zoom = (mask_xy &
             (np.abs(x_ec - x_zoom_c) < half_x + 0.5 * Dxz) &
             (np.abs(y_ec - y_zoom_c) < half_y + 0.5 * Dxz))
print(f'  Zoom patch: {mask_zoom.sum()} elements around '
      f'(x={x_zoom_c:.2f}, y={y_zoom_c:.2f})')

if mask_zoom.sum() == 0:
    print('  No elements in zoom region — skipping zoom figure.')
    raise SystemExit(0)

triang_zoom = build_triang(msh.x[mask_zoom, iz_mid, :, :],
                            msh.y[mask_zoom, iz_mid, :, :])
phi_zoom = phi_arr[mask_zoom, iz_mid, :, :].reshape(-1)
kap_zoom = last['kap_xy'].reshape(mask_xy.sum(), ly_, lx_)[
               mask_zoom[mask_xy], :, :].reshape(-1)

# Element boundary lines for zoom
x_zoom_lines = np.unique(np.round(msh.x[mask_zoom, iz_mid, :, 0].reshape(-1), 5))
y_zoom_lines = np.unique(np.round(msh.y[mask_zoom, iz_mid, 0, :].reshape(-1), 5))

# Quiver: element-local n̂ at all GLL points where |∇φ| > threshold
# Threshold: roughly 0.3 × peak gradient = 0.3 / (4ε) ≈ 0.075/ε
grad_threshold = 0.3 / (4.0 * EPS)

sub_phi = phi_arr[mask_zoom]
xs_z = msh.x[mask_zoom]; ys_z = msh.y[mask_zoom]; zs_z = msh.z[mask_zoom]
qx_all, qy_all, qu_all, qv_all = [], [], [], []
n_el_zoom = sub_phi.shape[0]
for e in range(n_el_zoom):
    x1d = xs_z[e, 0, 0, :]; y1d = ys_z[e, 0, :, 0]; z1d = zs_z[e, :, 0, 0]
    gz, gy, gx = np.gradient(sub_phi[e], z1d, y1d, x1d)
    gx2d = gx[iz_mid]; gy2d = gy[iz_mid]
    mag2d = np.sqrt(gx2d**2 + gy2d**2)
    mask_iface = mag2d > grad_threshold
    if not mask_iface.any():
        continue
    xx, yy = np.meshgrid(x1d, y1d, indexing='xy')
    mag_safe = np.where(mag2d > grad_threshold, mag2d, 1.0)
    nx2d = gx2d / mag_safe
    ny2d = gy2d / mag_safe
    qx_all.append(xx[mask_iface])
    qy_all.append(yy[mask_iface])
    qu_all.append(nx2d[mask_iface])
    qv_all.append(ny2d[mask_iface])

if qx_all:
    qx = np.concatenate(qx_all)
    qy = np.concatenate(qy_all)
    qu = np.concatenate(qu_all)
    qv = np.concatenate(qv_all)
    print(f'  Quiver: {len(qx)} interface-region GLL points '
          f'(|∇φ| > {grad_threshold:.1f})')
else:
    qx = qy = qu = qv = np.array([])
    print('  No interface GLL points found for quiver.')

arrow_scale = 0.6 * Dxz   # arrow length ≈ 0.6 element widths

fig2, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 5), dpi=args.dpi)
fig2.patch.set_facecolor('white')
fig2.subplots_adjust(left=0.06, right=0.96, top=0.88, bottom=0.10,
                     wspace=0.30)

for ax, vals, cmap, vmin, vmax, cb_label, title_str in [
    (ax_l, phi_zoom, 'RdBu_r', 0.0, 1.0, 'φ',
     f'φ + n̂ quiver   t = {last["t"]:.2f} TU'),
    (ax_r, kap_zoom, 'seismic', -kap_scale, kap_scale, 'κ',
     f'κ + n̂ quiver   κ_rms = {last["kappa_rms"]:.1f}  (ref {KAPPA_SPHERE:.1f})'),
]:
    ax.tricontourf(triang_zoom, vals, levels=LEVELS, cmap=cmap,
                   vmin=vmin, vmax=vmax, extend='both')
    sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig2.colorbar(sm, ax=ax, pad=0.02, aspect=14, fraction=0.04)
    cb.set_label(cb_label, fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # Thick element boundaries — kink locations
    for xv in x_zoom_lines:
        ax.axvline(xv, color='#333333', lw=0.6, alpha=0.7, zorder=1)
    for yv in y_zoom_lines:
        ax.axhline(yv, color='#333333', lw=0.6, alpha=0.7, zorder=1)

    # Interface contour
    ax.tricontour(triang_zoom, phi_zoom, levels=[0.5],
                  colors='limegreen', linewidths=1.8, zorder=3)

    # n̂ quiver — unit arrows, only at interface GLL points
    if len(qx) > 0:
        ax.quiver(qx, qy, qu, qv,
                  color='black', scale=1.0 / arrow_scale,
                  scale_units='xy', width=0.003,
                  headwidth=3, headlength=4, headaxislength=3,
                  zorder=4, alpha=0.85)

    ax.set_xlim(x_zoom_c - half_x - 0.5 * Dxz, x_zoom_c + half_x + 0.5 * Dxz)
    ax.set_ylim(y_zoom_c - half_y - 0.5 * Dxz, y_zoom_c + half_y + 0.5 * Dxz)
    ax.set_title(title_str, fontsize=10, fontweight='bold')
    ax.set_xlabel('x', fontsize=9)
    ax.set_ylabel('y', fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_aspect('equal')

fig2.suptitle(
    f'{run_label}   R={R}  ε={EPS}   zoom: top of drop  '
    f'(element size Δ≈{Dxz:.4f},  ε/Δ≈{EPS/Dxz:.2f},  ε/Δ_GLL≈{EPS/Dxz*7:.1f})',
    fontsize=10, fontweight='bold')

out_zoom = os.path.join(OUT_DIR, f'sigma0_normals_{args.run}.png')
plt.savefig(out_zoom, dpi=args.dpi, bbox_inches='tight', facecolor='white')
plt.close(fig2)
print(f'Saved: {out_zoom}')
print('Done.')
