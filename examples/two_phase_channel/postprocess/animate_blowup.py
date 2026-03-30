"""
animate_blowup.py — High-quality blow-up animations for supervisor presentation.

Produces three panels stacked vertically per frame:
  Top:    φ field (drop shape and interface integrity)
  Middle: κ (curvature, postprocess) — shows the artifact / blow-up cause
  Bottom: |u| (velocity magnitude — turbulent flow context)

Each frame is rendered as a fresh figure (colorbars always correct).
Frames are combined into a GIF via PIL.

Usage:
    python3 animate_blowup.py --run channel_test_laminar
    python3 animate_blowup.py --run channel_test_sigma0_diag --fps 3 --kappa-scale 35
    python3 animate_blowup.py --run channel_p3_sigma0 --mesh l2 --R 0.4 --fps 2
    python3 animate_blowup.py --run channel_l3_sigma0 --mesh l3 --R 0.4 --fps 2
"""

import argparse
import glob
import math
import os
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pandas as pd
from PIL import Image

from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.msh import Mesh as msh_c
from pysemtools.datatypes.field import Field as field_c

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--run', default='channel_test_laminar')
parser.add_argument('--fps', type=int, default=2)
parser.add_argument('--dpi', type=int, default=160)
parser.add_argument('--kappa-scale', type=float, default=None,
                    help='Symmetric κ colour limit (auto if not set)')
parser.add_argument('--no-gif', action='store_true')
parser.add_argument('--stride', type=int, default=1,
                    help='Animate every Nth snapshot (default: 1 = all)')
parser.add_argument('--mesh', choices=['p1', 'p2', 'l1', 'l2', 'l3', 'l4'], default=None,
                    help='Mesh preset (sets nz_elems for z-slice selection): '
                         'p1=81x18x27 (nz=27), p2/l1=108x18x36 (nz=36), '
                         'l2=144x24x48 (nz=48), l3=192x32x64 (nz=64), l4=288x48x96 (nz=96).')
parser.add_argument('--R', type=float, default=0.4,
                    help='Drop radius (default 0.4 for convergence series; old baseline uses 0.3)')
args = parser.parse_args()

RUN_DIR    = f'/lscratch/sieburgh/simulations/{args.run}'
OUT_DIR    = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
FRAMES_DIR = os.path.join(OUT_DIR, 'frames')
os.makedirs(FRAMES_DIR, exist_ok=True)
comm    = MPI.COMM_WORLD

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
LLX = 4.0 * math.pi
LLY = 2.0
LLZ = 4.0 / 3.0 * math.pi
Z_C = LLZ / 2.0
R   = args.R
KAPPA_SPHERE = 2.0 / R

_MESH_PRESETS = {'p1': 27, 'p2': 36, 'l1': 36, 'l2': 48, 'l3': 64, 'l4': 96}
_mesh_preset  = args.mesh if args.mesh else 'p1'
NZ_ELEMS_DEFAULT = _MESH_PRESETS[_mesh_preset]

# Spin-up run for mesh coordinates (two-phase restart files omit XYZ)
_SPINUP_RUNS = {
    'p1': 'channel_single_phase',
    'p2': 'channel_p2_single_phase', 'l1': 'channel_p2_single_phase',
    'l2': 'channel_p3_single_phase',
    'l3': 'channel_l3_single_phase',
    'l4': 'channel_l4_single_phase',
}
SIM_DIR = '/lscratch/sieburgh/simulations'

# ---------------------------------------------------------------------------
# Load ekin.csv
# ---------------------------------------------------------------------------
ekin_df = None
ekin_path = os.path.join(RUN_DIR, 'ekin.csv')
if os.path.exists(ekin_path):
    ekin_df = pd.read_csv(
        ekin_path, comment='#', header=None,
        names=['t', 'Ekin', 'enst', 'u_max',
               'kappa_max', 'kappa_min', 'kappa_rms',
               'Fst_max', 'phi_min', 'phi_max', 'extra'],
    )
    print(f'ekin.csv: {len(ekin_df)} rows')

def ekin_at(t):
    if ekin_df is None:
        return {}
    idx = (ekin_df['t'] - t).abs().idxmin()
    return ekin_df.iloc[idx].to_dict()

# ---------------------------------------------------------------------------
# Field files
# ---------------------------------------------------------------------------
field_files = sorted(glob.glob(os.path.join(RUN_DIR, 'field0.f[0-9]*')))
if not field_files:
    raise FileNotFoundError(f'No field0.f* in {RUN_DIR}')
if args.stride > 1:
    field_files = field_files[::args.stride]
print(f'Found {len(field_files)} field files in {RUN_DIR} (stride={args.stride})')

# ---------------------------------------------------------------------------
# Mesh  (two-phase restart files omit XYZ — fall back to spin-up f00000)
# ---------------------------------------------------------------------------
_spinup_f0 = os.path.join(SIM_DIR, _SPINUP_RUNS[_mesh_preset], 'field0.f00000')
_mesh_candidates = [os.path.join(RUN_DIR, 'field0.f00000'), _spinup_f0]
mesh_file = next((f for f in _mesh_candidates if os.path.exists(f)), field_files[0])
print(f'Reading mesh from: {mesh_file}')
xyz_data = preadnek(mesh_file, comm)
msh = msh_c(comm, data=xyz_data)

lz_ = msh.x.shape[1]
ly_ = msh.x.shape[2]
lx_ = msh.x.shape[3]
iz_mid = lz_ // 2
iy_mid = ly_ // 2

nz_elems = NZ_ELEMS_DEFAULT
dz_approx = LLZ / nz_elems
z_ec = msh.z[:, iz_mid, iy_mid, lx_ // 2]
mask_xy = np.abs(z_ec - Z_C) < 0.6 * dz_approx
print(f'  x-y slice (z≈{Z_C:.2f}): {mask_xy.sum()} elements')

# Triangulation for x-y slice
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

# Element boundary lines (for φ subplot)
x_lines = np.unique(np.round(msh.x[mask_xy, iz_mid, :, 0].reshape(-1), 5))
y_lines = np.unique(np.round(msh.y[mask_xy, iz_mid, 0, :].reshape(-1), 5))

# ---------------------------------------------------------------------------
# κ computation (element-local on x-y slice only — fast)
# ---------------------------------------------------------------------------
def compute_kappa_slice(phi_arr):
    sub = phi_arr[mask_xy]
    xs = msh.x[mask_xy];  ys = msh.y[mask_xy];  zs = msh.z[mask_xy]
    n_el = sub.shape[0]
    kappa_out = np.zeros((n_el, ly_, lx_))
    for e in range(n_el):
        x1d = xs[e, 0, 0, :];  y1d = ys[e, 0, :, 0];  z1d = zs[e, :, 0, 0]
        gz, gy, gx = np.gradient(sub[e], z1d, y1d, x1d)
        mag = np.sqrt(gx**2 + gy**2 + gz**2)
        mag = np.where(mag < 1e-3, 1e-3, mag)
        nx_ = gx / mag;  ny_ = gy / mag;  nz_ = gz / mag
        _, _, dnx = np.gradient(nx_, z1d, y1d, x1d)
        _, dny, _ = np.gradient(ny_, z1d, y1d, x1d)
        dnz, _, _ = np.gradient(nz_, z1d, y1d, x1d)
        kappa_out[e] = -(dnx + dny + dnz)[iz_mid]
    return kappa_out.reshape(-1)

# ---------------------------------------------------------------------------
# Load all snapshots
# ---------------------------------------------------------------------------
frames = []
print('Loading snapshots ...')
for fname in field_files:
    print(f'  {os.path.basename(fname)} ...', end=' ', flush=True)
    data = preadnek(fname, comm)
    fld  = field_c(comm, data=data)
    t    = fld.t
    phi_arr = fld.fields['scal'][0]
    vel_u   = fld.fields['vel'][0]
    vel_v   = fld.fields['vel'][1]
    vel_w   = fld.fields['vel'][2]
    phi_xy  = phi_arr[mask_xy, iz_mid, :, :].reshape(-1)
    kap_xy  = compute_kappa_slice(phi_arr)
    vel_xy  = np.sqrt(
        vel_u[mask_xy, iz_mid, :, :].reshape(-1)**2 +
        vel_v[mask_xy, iz_mid, :, :].reshape(-1)**2 +
        vel_w[mask_xy, iz_mid, :, :].reshape(-1)**2
    )
    ek      = ekin_at(t)
    frames.append({'t': t, 'phi_xy': phi_xy, 'kap_xy': kap_xy, 'vel_xy': vel_xy,
                   'kappa_rms_ek': ek.get('kappa_rms', np.nan),
                   'Fst_max':      ek.get('Fst_max', np.nan),
                   'phi_max':      float(phi_arr.reshape(-1).max())})
    print(f't={t:.3f}  κ_rms(Neko)={ek.get("kappa_rms", float("nan")):.1f}  '
          f'Fst_max={ek.get("Fst_max", float("nan")):.1f}')

# ---------------------------------------------------------------------------
# Colour scales — fixed across all frames
# ---------------------------------------------------------------------------
kap_absmax = max(np.percentile(np.abs(f['kap_xy']), 99) for f in frames)
vel_absmax = max(np.percentile(f['vel_xy'], 99) for f in frames)
kap_scale  = args.kappa_scale if args.kappa_scale else round(min(kap_absmax, 50.0))
vel_scale  = round(min(vel_absmax * 1.1, 10.0), 1)
print(f'\nColour limits:  κ ± {kap_scale:.0f},  |u| 0–{vel_scale:.1f}')

LEVELS = 100
run_name  = args.run
run_label = run_name.replace('channel_test_', '').replace('_', ' ')

# ---------------------------------------------------------------------------
# Render each frame as a standalone PNG
# ---------------------------------------------------------------------------
tmp_dir  = tempfile.mkdtemp()
png_paths = []

print('\nRendering frames ...')
for i, f in enumerate(frames):
    # Three panels stacked vertically: φ / κ / |u|
    # Wide landscape so the 6:1 domain fits without excess whitespace
    fig, (ax_phi, ax_kap, ax_vel) = plt.subplots(3, 1, figsize=(13, 7), dpi=args.dpi)
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(left=0.06, right=0.90, top=0.91, bottom=0.07, hspace=0.45)

    def _add_panel(ax, tcf_data, vals, cmap, vmin, vmax, title, cb_label, cb_ticks):
        # extend='both': out-of-range values render as saturated colours (not white).
        # This preserves visibility of simulation errors (e.g. phi << 0 or >> 1).
        # The colorbar is built from a ScalarMappable so it stays clean (no arrows).
        ax.tricontourf(triang_xy, vals, levels=LEVELS,
                       cmap=cmap, vmin=vmin, vmax=vmax, extend='both')
        ax.tricontour(triang_xy, f['phi_xy'], levels=[0.5],
                      colors='limegreen', linewidths=1.5)
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlabel('x', fontsize=8)
        ax.set_ylabel('y', fontsize=8)
        ax.tick_params(labelsize=7)
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, pad=0.01, aspect=12, fraction=0.025)
        cb.set_label(cb_label, fontsize=9)
        cb.set_ticks(cb_ticks)
        cb.ax.tick_params(labelsize=8)

    # ---- φ panel ----
    for xv in x_lines:
        ax_phi.axvline(xv, color='gray', lw=0.2, alpha=0.35, zorder=0)
    for yv in y_lines:
        ax_phi.axhline(yv, color='gray', lw=0.2, alpha=0.35, zorder=0)
    _add_panel(ax_phi, None, f['phi_xy'], 'RdBu_r', 0.0, 1.0,
               'Phase field  φ  (green = interface φ = 0.5)',
               'φ', [0, 0.5, 1.0])

    # ---- κ panel ----
    _add_panel(ax_kap, None, f['kap_xy'], 'seismic', -kap_scale, kap_scale,
               f'Curvature  κ  (sphere ref = {KAPPA_SPHERE:.1f},  colour limit ±{kap_scale:.0f})',
               'κ', [-kap_scale, 0, kap_scale])

    # ---- |u| panel ----
    _add_panel(ax_vel, None, f['vel_xy'], 'viridis', 0.0, vel_scale,
               f'Velocity magnitude  |u|  (colour limit 0–{vel_scale:.1f})',
               '|u|', [0, round(vel_scale / 2, 1), vel_scale])

    # ---- title / annotation ----
    kappa_str = f'κ_rms (Neko) = {f["kappa_rms_ek"]:.1f}' if not math.isnan(f['kappa_rms_ek']) else ''
    fst_str   = (f'   |F_ST|_max = {f["Fst_max"]:.1f}' if not math.isnan(f['Fst_max']) and f['Fst_max'] > 0
                 else '')
    fig.suptitle(
        f'{run_label}   t = {f["t"]:.4f} TU     frame {i+1}/{len(frames)}\n'
        f'{kappa_str}{fst_str}',
        fontsize=11, y=0.98, fontweight='bold'
    )

    out_png = os.path.join(tmp_dir, f'frame_{i:03d}.png')
    plt.savefig(out_png, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    png_paths.append(out_png)
    print(f'  Frame {i+1:02d}/{len(frames)}  t={f["t"]:.4f}')

# Copy frames to figures/frames/ with run name prefix
import shutil
for src in png_paths:
    basename = os.path.basename(src).replace('frame_', f'blowup_{run_name}_frame')
    shutil.copy(src, os.path.join(FRAMES_DIR, basename))

# ---------------------------------------------------------------------------
# Assemble GIF from PNGs
# ---------------------------------------------------------------------------
if not args.no_gif:
    out_gif = os.path.join(OUT_DIR, f'blowup_{run_name}.gif')
    print(f'\nAssembling GIF → {out_gif}')
    duration_ms = int(1000 / args.fps)
    imgs = [Image.open(p) for p in png_paths]
    imgs[0].save(
        out_gif,
        save_all=True,
        append_images=imgs[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
    )
    print(f'Done: {out_gif}  ({len(imgs)} frames, {args.fps} fps)')

import shutil
shutil.rmtree(tmp_dir, ignore_errors=True)
print('All done.')
