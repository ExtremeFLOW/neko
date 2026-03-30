"""
animate_three_meshes.py — Mesh convergence animation for σ=0 CDI runs.

Three rows stacked vertically per frame, each showing the x-y midplane
φ field for one mesh level:
  Row 0: L1  ε=0.053  108×18×36
  Row 1: L2  ε=0.040  144×24×48
  Row 2: L3  ε=0.030  192×32×64

Element boundary lines are drawn in each row — the different grid densities
make mesh refinement immediately visible. The φ=0.5 interface contour is
overlaid in each panel. κ_rms (from ekin.csv) is annotated in the title
of each row.

Frames are matched by nearest simulation time across the three runs, then
assembled into a GIF.

Data size note: L3 field files are 3.8 GB each. Ten sparse snapshots
(≈0.5 TU interval) take ≈60 s to read on egidius local XFS.

Usage:
    cd examples/two_phase_channel
    source ../../setup-env-channel.sh --egidius
    python3 postprocess/animate_three_meshes.py
    python3 postprocess/animate_three_meshes.py --fps 3 --stride 2
    python3 postprocess/animate_three_meshes.py --no-gif   # frames only
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
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from PIL import Image

from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.msh import Mesh as msh_c
from pysemtools.datatypes.field import Field as field_c

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--fps',    type=int,   default=2)
parser.add_argument('--dpi',    type=int,   default=160)
parser.add_argument('--stride', type=int,   default=1,
                    help='Use every Nth snapshot (default: 1 = all synced files)')
parser.add_argument('--no-gif', action='store_true',
                    help='Write frames but skip GIF assembly')
parser.add_argument('--kappa-scale', type=float, default=None,
                    help='Manual κ colour limit (auto if unset)')
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Mesh level configurations
# ---------------------------------------------------------------------------
SIM_DIR = '/lscratch/sieburgh/simulations'
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)
FRAMES_DIR = os.path.join(OUT_DIR, 'frames')
os.makedirs(FRAMES_DIR, exist_ok=True)

comm = MPI.COMM_WORLD

CONFIGS = [
    dict(run='channel_p2_sigma0_eps053', label='L1', eps=0.053, nz_elems=36,
         mesh_str='108×18×36', spinup='channel_p2_single_phase'),
    dict(run='channel_p3_sigma0',        label='L2', eps=0.040, nz_elems=48,
         mesh_str='144×24×48', spinup='channel_p3_single_phase'),
    dict(run='channel_l3_sigma0',        label='L3', eps=0.030, nz_elems=64,
         mesh_str='192×32×64', spinup='channel_l3_single_phase'),
]

R            = 0.4
KAPPA_SPHERE = 2.0 / R
LLX  = 4.0 * math.pi
LLZ  = 4.0 / 3.0 * math.pi
Z_C  = LLZ / 2.0
LEVELS = 80

# ---------------------------------------------------------------------------
# Load all mesh levels
# ---------------------------------------------------------------------------
def _build_triang(xs, ys):
    """Build matplotlib Triangulation from GLL element arrays."""
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


def load_level(cfg):
    """Load all snapshots for one mesh level. Returns a dataset dict."""
    run_dir = os.path.join(SIM_DIR, cfg['run'])
    print(f'\n=== Loading {cfg["label"]}  ({cfg["run"]}) ===')

    # Field files
    field_files = sorted(glob.glob(os.path.join(run_dir, 'field0.f[0-9]*')))
    if not field_files:
        raise FileNotFoundError(f'No field0.f* in {run_dir}')
    if args.stride > 1:
        field_files = field_files[::args.stride]
    print(f'  {len(field_files)} snapshots (stride={args.stride})')

    # ekin.csv
    ekin_df = None
    ekin_path = os.path.join(run_dir, 'ekin.csv')
    if os.path.exists(ekin_path):
        ekin_df = pd.read_csv(
            ekin_path, comment='#', header=None,
            names=['t', 'Ekin', 'enst', 'u_max',
                   'kappa_max', 'kappa_min', 'kappa_rms',
                   'Fst_max', 'phi_min', 'phi_max', 'extra'],
        )

    def kappa_rms_at(t):
        if ekin_df is None:
            return float('nan')
        idx = (ekin_df['t'] - t).abs().idxmin()
        return float(ekin_df.loc[idx, 'kappa_rms'])

    # Mesh: two-phase restart files don't include XYZ (var[0]=0).
    # Use spin-up field0.f00000 which has geometry embedded.
    _spinup_f0 = os.path.join(SIM_DIR, cfg['spinup'], 'field0.f00000')
    _candidates = [os.path.join(run_dir, 'field0.f00000'), _spinup_f0]
    mesh_file = next((f for f in _candidates if os.path.exists(f)), field_files[0])
    print(f'  Reading mesh from: {os.path.relpath(mesh_file, SIM_DIR)}')
    xyz_data = preadnek(mesh_file, comm)
    msh = msh_c(comm, data=xyz_data)

    lz_ = msh.x.shape[1]
    ly_ = msh.x.shape[2]
    lx_ = msh.x.shape[3]
    iz_mid = lz_ // 2

    dz_approx = LLZ / cfg['nz_elems']
    z_ec = msh.z[:, iz_mid, ly_ // 2, lx_ // 2]
    mask_xy = np.abs(z_ec - Z_C) < 0.6 * dz_approx
    print(f'  x-y slice: {mask_xy.sum()} elements')

    triang_xy = _build_triang(msh.x[mask_xy, iz_mid, :, :],
                               msh.y[mask_xy, iz_mid, :, :])

    # Element boundary coordinates
    x_lines = np.unique(np.round(msh.x[mask_xy, iz_mid, :, 0].reshape(-1), 5))
    y_lines = np.unique(np.round(msh.y[mask_xy, iz_mid, 0, :].reshape(-1), 5))

    # Snapshots
    frames = []
    for fname in field_files:
        print(f'  {os.path.basename(fname)} ...', end=' ', flush=True)
        data = preadnek(fname, comm)
        fld  = field_c(comm, data=data)
        t    = fld.t
        phi_xy = fld.fields['scal'][0][mask_xy, iz_mid, :, :].reshape(-1)
        kr = kappa_rms_at(t)
        frames.append({'t': t, 'phi_xy': phi_xy, 'kappa_rms': kr})
        print(f't={t:.3f}  κ_rms={kr:.1f}')

    return dict(cfg=cfg, triang=triang_xy, x_lines=x_lines, y_lines=y_lines,
                frames=frames, lx=lx_, ly=ly_)


datasets = [load_level(cfg) for cfg in CONFIGS]

# ---------------------------------------------------------------------------
# Match frames across mesh levels by nearest simulation time
# ---------------------------------------------------------------------------
ref_times = [f['t'] for f in datasets[0]['frames']]
matched = []
for t_ref in ref_times:
    row = []
    for ds in datasets:
        best = min(ds['frames'], key=lambda f: abs(f['t'] - t_ref))
        row.append(best)
    matched.append(row)

print(f'\n{len(matched)} matched frames  '
      f'(t = {ref_times[0]:.3f} → {ref_times[-1]:.3f} TU)')

# ---------------------------------------------------------------------------
# Render frames
# ---------------------------------------------------------------------------
# Element line style: slightly lighter alpha at finer meshes to avoid clutter
LINE_STYLES = [
    dict(color='#888888', lw=0.22, alpha=0.40),  # L1 — fewest lines
    dict(color='#888888', lw=0.18, alpha=0.35),  # L2
    dict(color='#888888', lw=0.15, alpha=0.28),  # L3 — most lines
]

tmp_dir  = tempfile.mkdtemp()
png_paths = []

print('\nRendering frames ...')
for fi, frame_row in enumerate(matched):
    t_ref = frame_row[0]['t']

    # Figure: 3 rows, one per mesh level
    # Wide format — 4π:2 ≈ 6.3:1 domain per row, thin panels
    fig, axes = plt.subplots(3, 1, figsize=(14, 7.5), dpi=args.dpi)
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(left=0.04, right=0.88, top=0.91, bottom=0.06,
                        hspace=0.50)

    for row_idx, (ax, frame, ds, ls) in enumerate(
            zip(axes, frame_row, datasets, LINE_STYLES)):
        cfg = ds['cfg']

        # φ filled contour
        ax.tricontourf(ds['triang'], frame['phi_xy'],
                       levels=LEVELS, cmap='RdBu_r', vmin=0.0, vmax=1.0,
                       extend='both')

        # Interface contour
        ax.tricontour(ds['triang'], frame['phi_xy'], levels=[0.5],
                      colors='limegreen', linewidths=1.4)

        # Element boundary lines
        for xv in ds['x_lines']:
            ax.axvline(xv, **ls, zorder=0)
        for yv in ds['y_lines']:
            ax.axhline(yv, **ls, zorder=0)

        # Axis formatting
        ax.set_xlim(0, LLX)
        ax.set_ylim(-1.0, 1.0)
        ax.set_ylabel('y', fontsize=8)
        ax.tick_params(labelsize=7)
        if row_idx == 2:
            ax.set_xlabel('x', fontsize=8)

        kr = frame['kappa_rms']
        kr_str = f'  κ_rms={kr:.1f}  (ref 2/R={KAPPA_SPHERE:.1f})' \
                 if not math.isnan(kr) else ''
        ax.set_title(
            f'{cfg["label"]}  ε={cfg["eps"]}  {cfg["mesh_str"]}{kr_str}',
            fontsize=9, fontweight='bold', loc='left', pad=3)

    # Shared φ colorbar on the right
    sm = ScalarMappable(cmap='RdBu_r', norm=Normalize(vmin=0.0, vmax=1.0))
    sm.set_array([])
    cbar_ax = fig.add_axes([0.90, 0.10, 0.015, 0.78])
    cb = fig.colorbar(sm, cax=cbar_ax)
    cb.set_label('φ', fontsize=10)
    cb.set_ticks([0.0, 0.5, 1.0])
    cb.ax.tick_params(labelsize=8)

    fig.suptitle(
        f'σ=0 CDI diagnostic — mesh convergence\n'
        f't = {t_ref:.3f} TU     frame {fi+1}/{len(matched)}',
        fontsize=11, y=0.98, fontweight='bold')

    out_png = os.path.join(tmp_dir, f'frame_{fi:03d}.png')
    plt.savefig(out_png, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    png_paths.append(out_png)
    print(f'  Frame {fi+1:02d}/{len(matched)}  t={t_ref:.3f}')

# Copy frames to figures/frames/
import shutil
for src in png_paths:
    dst = os.path.join(FRAMES_DIR,
                       os.path.basename(src).replace('frame_', 'three_meshes_frame'))
    shutil.copy(src, dst)

# ---------------------------------------------------------------------------
# Assemble GIF
# ---------------------------------------------------------------------------
if not args.no_gif:
    out_gif = os.path.join(OUT_DIR, 'three_meshes_sigma0.gif')
    print(f'\nAssembling GIF → {out_gif}')
    duration_ms = int(1000 / args.fps)
    imgs = [Image.open(p) for p in png_paths]
    imgs[0].save(out_gif, save_all=True, append_images=imgs[1:],
                 duration=duration_ms, loop=0, optimize=False)
    print(f'Done: {out_gif}  ({len(imgs)} frames, {args.fps} fps)')

shutil.rmtree(tmp_dir, ignore_errors=True)
print('All done.')
