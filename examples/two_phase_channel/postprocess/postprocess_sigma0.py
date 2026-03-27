"""
postprocess_sigma0.py — Clean diagnostic figures for σ=0 CDI runs.

Replaces analyze_sigma0_normals.py with a faster, better-looking script.

Produces:
  1. Time-series figure  (reads ekin.csv — always fast)
        κ_rms(t)  with 2/R reference line
        φ_max(t) and φ_min(t)
        E_kin(t)  — turbulence stability check

  2. Field snapshots     (optional — loads 3 field files: first / mid / last)
        x-y slice: φ (top row) and κ (bottom row)
        Same rendering style as animate_blowup.py

Usage:
    python3 postprocess_sigma0.py --run channel_test_sigma0_long
    python3 postprocess_sigma0.py --run channel_test_sigma0_long --no-snapshots
    python3 postprocess_sigma0.py --run channel_p2_sigma0_eps053 \\
        --R 0.4 --eps 0.053 --mesh p2
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
                    help='Interface width ε (default: 0.07 for legacy 81x18x27)')
parser.add_argument('--mesh', choices=['p1', 'p2', 'l2', 'l3'], default='p1',
                    help='Mesh preset for z-slice selection '
                         '(p1=81x18x27 nz=27, p2=108x18x36 nz=36, '
                         'l2=144x24x48 nz=48, l3=192x32x64 nz=64)')
parser.add_argument('--kappa-scale', type=float, default=None,
                    help='Symmetric κ colour limit for snapshots (auto if unset)')
parser.add_argument('--no-snapshots', action='store_true',
                    help='Skip field file loading; produce only the time-series figure')
parser.add_argument('--dpi', type=int, default=160)
args = parser.parse_args()

RUN_DIR = f'/lscratch/sieburgh/simulations/{args.run}'
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
comm    = MPI.COMM_WORLD

R            = args.R
EPS          = args.eps
KAPPA_SPHERE = 2.0 / R

LLZ   = 4.0 / 3.0 * math.pi
Z_C   = LLZ / 2.0
_NZ   = {'p1': 27, 'p2': 36, 'l2': 48, 'l3': 64}
NZ_ELEMS = _NZ[args.mesh]

run_label = args.run.replace('channel_test_', '').replace('channel_', '').replace('_', ' ')

# ---------------------------------------------------------------------------
# 1. Time-series from ekin.csv (fast — no field files needed)
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

# Panel 1: κ_rms
ax1.plot(df.t, df.kappa_rms, color='steelblue', lw=1.5)
ax1.axhline(KAPPA_SPHERE, color='crimson', ls='--', lw=1.2,
            label=f'2/R = {KAPPA_SPHERE:.2f}  (sphere)')
ax1.set_ylabel('κ_rms', fontsize=11)
ax1.set_xlabel('t  [TU]', fontsize=10)
ax1.set_title('Curvature rms', fontsize=11, fontweight='bold')
ax1.legend(fontsize=9)
ax1.tick_params(labelsize=9)
ax1.set_xlim(df.t.min(), df.t.max())

# Panel 2: φ_max and φ_min
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

# Panel 3: E_kin
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
# 2. Field snapshots: first / mid / last  (optional)
# ---------------------------------------------------------------------------
if args.no_snapshots:
    print('Skipping field snapshots (--no-snapshots).')
    raise SystemExit(0)

field_files = sorted(glob.glob(os.path.join(RUN_DIR, 'field0.f[0-9]*')))
if not field_files:
    print('No field files found — skipping snapshots.')
    raise SystemExit(0)

n = len(field_files)
idx_first = 0
idx_mid   = n // 2
idx_last  = n - 1
snap_files = [field_files[i] for i in sorted({idx_first, idx_mid, idx_last})]
print(f'\nLoading {len(snap_files)} snapshots for field plots ...')

# --- Mesh (from first field file) ---
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

# Triangulation
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

# κ computation (element-local, x-y slice only)
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

# Load snapshots
snaps = []
for fpath in snap_files:
    print(f'  {os.path.basename(fpath)} ...', end=' ', flush=True)
    data = preadnek(fpath, comm)
    fld  = field_c(comm, data=data)
    t    = fld.t
    phi_arr = fld.fields['scal'][0]
    phi_xy  = phi_arr[mask_xy, iz_mid, :, :].reshape(-1)
    kap_xy  = compute_kappa_slice(phi_arr)
    row = df.iloc[(df.t - t).abs().argsort().iloc[0]]
    print(f't={t:.3f}  κ_rms={row.kappa_rms:.1f}  φ_max={row.phi_max:.4f}')
    snaps.append({'t': t, 'phi_xy': phi_xy, 'kap_xy': kap_xy,
                  'kappa_rms': row.kappa_rms})

# Colour scale: fixed across all snapshots
kap_absmax = max(np.percentile(np.abs(s['kap_xy']), 99) for s in snaps)
kap_scale  = args.kappa_scale if args.kappa_scale else round(min(kap_absmax, 50.0))
print(f'\nκ colour limit: ±{kap_scale:.0f}')

# --- Render snapshots figure (2 rows × N_snaps columns) ---
LEVELS = 100
N      = len(snaps)
fig, axes = plt.subplots(2, N, figsize=(5 * N, 6), dpi=args.dpi)
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.05, right=0.93, top=0.90, bottom=0.07,
                    hspace=0.35, wspace=0.25)

if N == 1:
    axes = axes.reshape(2, 1)

def _add_panel(ax, triang, vals, cmap, vmin, vmax, title, cb_label, cb_ticks):
    ax.tricontourf(triang, vals, levels=LEVELS, cmap=cmap,
                   vmin=vmin, vmax=vmax, extend='both')
    sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02, aspect=14, fraction=0.03)
    cb.set_label(cb_label, fontsize=8)
    cb.set_ticks(cb_ticks)
    cb.ax.tick_params(labelsize=7)
    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.set_xlabel('x', fontsize=8)
    ax.set_ylabel('y', fontsize=8)
    ax.tick_params(labelsize=7)

for col, s in enumerate(snaps):
    # φ panel — with element boundaries
    ax_phi = axes[0, col]
    for xv in x_lines:
        ax_phi.axvline(xv, color='gray', lw=0.2, alpha=0.35, zorder=0)
    for yv in y_lines:
        ax_phi.axhline(yv, color='gray', lw=0.2, alpha=0.35, zorder=0)
    _add_panel(ax_phi, triang_xy, s['phi_xy'], 'RdBu_r', 0.0, 1.0,
               f'φ   t = {s["t"]:.3f} TU', 'φ', [0, 0.5, 1.0])
    ax_phi.tricontour(triang_xy, s['phi_xy'], levels=[0.5],
                      colors='limegreen', linewidths=1.5)

    # κ panel
    ax_kap = axes[1, col]
    _add_panel(ax_kap, triang_xy, s['kap_xy'], 'seismic',
               -kap_scale, kap_scale,
               f'κ   κ_rms = {s["kappa_rms"]:.1f}  (ref {KAPPA_SPHERE:.1f})',
               'κ', [-kap_scale, 0, kap_scale])
    ax_kap.tricontour(triang_xy, s['phi_xy'], levels=[0.5],
                      colors='limegreen', linewidths=1.5)

fig.suptitle(f'{run_label}   R={R}  ε={EPS}   field snapshots  (first / mid / last)',
             fontsize=11, fontweight='bold')

out_snaps = os.path.join(OUT_DIR, f'sigma0_snapshots_{args.run}.png')
plt.savefig(out_snaps, dpi=args.dpi, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out_snaps}')
print('Done.')
