"""
Animation of the two-phase turbulent channel flow.

Requires pysemtools (pip install pysemtools).

Shows two views of the phase field and velocity magnitude:
  Left:  x-y side-view slice at z ≈ z_c (streamwise × wall-normal)
  Right: x-z top-view  slice at y ≈ 0   (streamwise × spanwise)

Both views show the phi=0.5 drop interface as a contour line.

Usage:
    python3 animate_two_phase_channel.py
    python3 animate_two_phase_channel.py --run channel_test_v2 --fps 4

Output: animate_<run_name>.gif  (saved next to this script)
"""

import argparse
import glob
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np

from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.msh import Mesh as msh_c
from pysemtools.datatypes.field import Field as field_c

# -------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--run', default='channel_test_v4',
                    help='Run directory name under /lscratch/sieburgh/simulations/')
parser.add_argument('--fps', type=int, default=3,
                    help='Frames per second in output GIF')
parser.add_argument('--dpi', type=int, default=120)
parser.add_argument('--nz_elems', type=int, default=27,
                    help='Number of mesh elements in z (18 for 18³ runs, 27 for 81×18×27)')
args = parser.parse_args()

RUN_DIR  = f'/lscratch/sieburgh/simulations/{args.run}'
OUT_DIR  = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
OUT_FILE = os.path.join(OUT_DIR, f'animate_{args.run}.gif')

comm = MPI.COMM_WORLD

# -------------------------------------------------------------------
# Domain constants (from turb_channel_two_phase.f90)
# -------------------------------------------------------------------
LLX = 4.0 * math.pi          # streamwise length ≈ 12.57
LLZ = 4.0 / 3.0 * math.pi   # spanwise  length  ≈  4.19
X_C = LLX / 2.0              # drop centre x  ≈ 6.28
Z_C = LLZ / 2.0              # drop centre z  ≈ 2.09

# -------------------------------------------------------------------
# Discover field files
# -------------------------------------------------------------------
field_files = sorted(glob.glob(os.path.join(RUN_DIR, 'field0.f[0-9]*')))
if not field_files:
    raise FileNotFoundError(f'No field0.f* files in {RUN_DIR}')
print(f'Found {len(field_files)} field files in {RUN_DIR}')

# -------------------------------------------------------------------
# Read mesh from the first file and build slice triangulations
# -------------------------------------------------------------------
print('Reading mesh ...')
xyz_info = preadnek(field_files[0], comm)
msh = msh_c(comm, data=xyz_info)

n_elem = msh.x.shape[0]
lz     = msh.x.shape[1]
ly     = msh.x.shape[2]
lx     = msh.x.shape[3]
iz_mid = lz // 2   # GLL z-index used for x-y slice (index 3 for lz=6)
iy_mid = ly // 2   # GLL y-index used for x-z slice

print(f'  n_elem={n_elem}, lx=ly=lz={lx}  (iz_mid={iz_mid}, iy_mid={iy_mid})')

# z-coordinate at the "centre" GLL point of each element
z_elem_centers = msh.z[:, iz_mid, iy_mid, lx // 2]
y_elem_centers = msh.y[:, iz_mid, iy_mid, lx // 2]

# Approximate element widths (uniform in x and z)
nz_elems  = args.nz_elems
dz_approx = LLZ / nz_elems  # elements in z
dy_approx = 2.0  / 18.0  # 18 elements in y (max, near centre)

# Masks for the two slices
mask_xy = np.abs(z_elem_centers - Z_C) < 0.6 * dz_approx   # x-y at z≈Z_C
mask_xz = np.abs(y_elem_centers)        < 0.6 * dy_approx   # x-z at y≈0

print(f'  x-y slice (z≈{Z_C:.2f}): {mask_xy.sum()} elements')
print(f'  x-z slice (y≈0):        {mask_xz.sum()} elements')


def build_triang(xs, ys):
    """Build a matplotlib Triangulation from a set of (x,y) GLL points."""
    n_el, nj, ni = xs.shape
    pts_x, pts_y, tris = [], [], []
    offset = 0
    for e in range(n_el):
        pts_x.append(xs[e].flatten())
        pts_y.append(ys[e].flatten())
        for j in range(nj - 1):
            for i in range(ni - 1):
                p0 = offset + j * ni + i
                p1 = p0 + 1
                p2 = offset + (j + 1) * ni + i
                p3 = p2 + 1
                tris.append([p0, p1, p3])
                tris.append([p0, p3, p2])
        offset += nj * ni
    return tri.Triangulation(
        np.concatenate(pts_x),
        np.concatenate(pts_y),
        np.array(tris),
    )


# x-y triangulation (plot coords: x horizontal, y vertical)
triang_xy = build_triang(
    msh.x[mask_xy, iz_mid, :, :],
    msh.y[mask_xy, iz_mid, :, :],
)
# x-z triangulation (plot coords: x horizontal, z vertical)
triang_xz = build_triang(
    msh.x[mask_xz, :, iy_mid, :],
    msh.z[mask_xz, :, iy_mid, :],
)

# -------------------------------------------------------------------
# Load all snapshots
# -------------------------------------------------------------------
print('Loading snapshots ...')
snaps = []
for fname in field_files:
    data = preadnek(fname, comm)
    fld  = field_c(comm, data=data)

    phi_all = fld.fields['scal'][0]
    u_all   = fld.fields['vel'][0]
    v_all   = fld.fields['vel'][1]
    w_all   = fld.fields['vel'][2]
    umag_all = np.sqrt(u_all**2 + v_all**2 + w_all**2)

    # x-y slice values
    phi_xy  = phi_all [mask_xy, iz_mid, :, :].reshape(-1)
    umag_xy = umag_all[mask_xy, iz_mid, :, :].reshape(-1)

    # x-z slice values
    phi_xz  = phi_all [mask_xz, :, iy_mid, :].reshape(-1)
    umag_xz = umag_all[mask_xz, :, iy_mid, :].reshape(-1)

    snaps.append({
        't':      fld.t,
        'phi_xy':  phi_xy,
        'umag_xy': umag_xy,
        'phi_xz':  phi_xz,
        'umag_xz': umag_xz,
    })
    print(f'  t={fld.t:.3f}  phi_max={phi_all.max():.3f}')

# Global colour scales (fixed across all frames)
vel_max  = max(s['umag_xy'].max() for s in snaps)
phi_min_ = min(s['phi_xy'].min()  for s in snaps)

print(f'\nGlobal |u|_max = {vel_max:.4f}')

# -------------------------------------------------------------------
# Build figure (2 rows × 2 columns: phi / |u| × xy / xz)
# -------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 6), dpi=args.dpi)
fig.subplots_adjust(right=0.88, wspace=0.08, hspace=0.25)

ax_phi_xy, ax_phi_xz = axes[0]
ax_vel_xy, ax_vel_xz = axes[1]

for ax in axes.flat:
    ax.set_aspect('equal')

ax_phi_xy.set(title='φ  |  x–y (z=z_c)', xlabel='x', ylabel='y')
ax_phi_xz.set(title='φ  |  x–z (y=0)',   xlabel='x', ylabel='z')
ax_vel_xy.set(title='|u|  |  x–y',        xlabel='x', ylabel='y')
ax_vel_xz.set(title='|u|  |  x–z',        xlabel='x', ylabel='z')

# Fixed colorbars
cax_phi = fig.add_axes([0.89, 0.54, 0.015, 0.38])
cax_vel = fig.add_axes([0.89, 0.08, 0.015, 0.38])
sm_phi = ScalarMappable(cmap='RdBu_r',  norm=Normalize(vmin=0, vmax=1))
sm_vel = ScalarMappable(cmap='viridis', norm=Normalize(vmin=0, vmax=vel_max))
sm_phi.set_array([]); sm_vel.set_array([])
fig.colorbar(sm_phi, cax=cax_phi, label='φ')
fig.colorbar(sm_vel, cax=cax_vel, label='|u|')

title_obj = fig.suptitle('', fontsize=11)


def _update(frame):
    s = snaps[frame]
    for ax in axes.flat:
        for coll in ax.collections:
            coll.remove()

    for ax, trg, vals, cmap, vmax in [
        (ax_phi_xy, triang_xy, s['phi_xy'],  'RdBu_r',  1.0),
        (ax_phi_xz, triang_xz, s['phi_xz'],  'RdBu_r',  1.0),
        (ax_vel_xy, triang_xy, s['umag_xy'], 'viridis', vel_max),
        (ax_vel_xz, triang_xz, s['umag_xz'], 'viridis', vel_max),
    ]:
        ax.tricontourf(trg, vals, levels=80, cmap=cmap, vmin=0, vmax=vmax, extend='both')

    # phi=0.5 interface contour
    for ax, trg, phi_vals, colour in [
        (ax_phi_xy, triang_xy, s['phi_xy'],  'k'),
        (ax_phi_xz, triang_xz, s['phi_xz'],  'k'),
        (ax_vel_xy, triang_xy, s['phi_xy'],  'w'),
        (ax_vel_xz, triang_xz, s['phi_xz'],  'w'),
    ]:
        try:
            ax.tricontour(trg, phi_vals, levels=[0.5], colors=colour, linewidths=1.2)
        except Exception:
            pass   # contour may not exist if drop has dissolved

    phi_now = snaps[frame]['phi_xy'].max()
    title_obj.set_text(f't = {s["t"]:.2f}   |   φ_max = {phi_now:.3f}')


# -------------------------------------------------------------------
# Render animation
# -------------------------------------------------------------------
print(f'\nRendering {len(snaps)}-frame animation ...')
anim = FuncAnimation(fig, _update, frames=len(snaps),
                     interval=int(1000 / args.fps), repeat=True)
anim.save(OUT_FILE, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
plt.close(fig)
print(f'Saved: {OUT_FILE}')
