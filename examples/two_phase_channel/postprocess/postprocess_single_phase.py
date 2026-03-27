#!/usr/bin/env python3
"""
Single-phase turbulent channel — postprocessing.

Produces two figures:
  1. ekin_single_phase.png  — u_max and E_kin vs t (spin-up / turbulence indicator)
  2. meanprofile_single_phase.png — time-averaged mean velocity profile vs DNS / Reichardt

Usage (serial — no mpirun needed):
    source ../../setup-env-channel.sh --egidius
    python3 postprocess_single_phase.py

The mean profile is computed by time-averaging the last NSNAP field snapshots
(default: f00020–f00025, i.e. t=20–25 TU) and averaging over the homogeneous
streamwise (x) and spanwise (z) directions.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ── paths ─────────────────────────────────────────────────────────────────────
RUN_DIR = '/lscratch/sieburgh/simulations/channel_single_phase'
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')

Re_tau = 180.0
Re_b   = 2800.0
u_tau  = Re_tau / Re_b   # in bulk-velocity units (U_b=1)

# ── helper: Reichardt mean profile ────────────────────────────────────────────
def reichardt(yp, k=0.41, C=5.17):
    """u+ as a function of y+ using the Reichardt formula."""
    return (1/k * np.log(1 + k*yp)
            + (C - np.log(k)/k) * (1 - np.exp(-yp/11) - yp/11*np.exp(-yp/3)))

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — ekin.csv time series
# ─────────────────────────────────────────────────────────────────────────────
ekin_path = os.path.join(RUN_DIR, 'ekin.csv')
data = np.loadtxt(ekin_path, delimiter=',', comments='#')
t      = data[:, 0]
Ekin   = data[:, 1]
u_max  = data[:, 2]

# Reichardt centerline: u_cl = u^+_cl * u_tau
yp_cl   = Re_tau   # centerline y+ = h * Re_tau / h = Re_tau
u_cl    = reichardt(yp_cl) * u_tau  # in U_b units

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# u_max vs t
ax = axes[0]
ax.plot(t, u_max, color='#1f4e79', lw=1.3, label='simulation')
ax.axhline(1.5, color='gray', lw=1.0, ls='--', label='Poiseuille ($U_{cl}/U_b = 3/2$)')
ax.axhline(u_cl, color='#cc2222', lw=1.0, ls='--',
           label=f'Reichardt $U_{{cl}}^+·u_\\tau = {u_cl:.3f}$')
ax.axvline(20.0, color='green', lw=1.0, ls=':', alpha=0.8,
           label='checkpoint t=20 (fluid00004.chkp)')
ax.set(xlabel='$t$ (TU)', ylabel='$u_{max}$',
       title='Max velocity — spin-up to turbulence')
ax.legend(fontsize=8)

# E_kin vs t
ax = axes[1]
ax.plot(t, Ekin, color='#1f4e79', lw=1.3)
ax.axvline(20.0, color='green', lw=1.0, ls=':', alpha=0.8, label='t=20 (checkpoint)')
ax.set(xlabel='$t$ (TU)', ylabel='$E_{kin}$',
       title='Kinetic energy')
ax.legend(fontsize=8)

fig.suptitle(f'Single-phase channel  ($Re_b={Re_b:.0f}$, $Re_\\tau={Re_tau:.0f}$)',
             fontsize=12)
plt.tight_layout()
out1 = os.path.join(OUT_DIR, 'ekin_single_phase.png')
plt.savefig(out1, dpi=150)
plt.close()
print(f'Saved: {out1}')

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — mean velocity profile from field snapshots
# ─────────────────────────────────────────────────────────────────────────────
try:
    from mpi4py import MPI
    from pysemtools.io.ppymech.neksuite import preadnek
except ImportError as e:
    print(f'pysemtools not available ({e}); skipping mean profile figure.')
    raise SystemExit(0)

comm = MPI.COMM_WORLD

# ── read mesh coordinates from initial file (only f00000 has XYZ) ────────────
mesh_file = os.path.join(RUN_DIR, 'field0.f00000')
print(f'Reading mesh from {os.path.basename(mesh_file)} ...')
xyz_data = preadnek(mesh_file, comm)
# Coordinates: elem.pos shape (3, lz, ly, lx); pos[1] = wall-normal y
y_all = np.concatenate([e.pos[1].flatten() for e in xyz_data.elem])
print(f'  {len(xyz_data.elem)} elements, y ∈ [{y_all.min():.3f}, {y_all.max():.3f}]')

# ── indices of snapshots to average (t ≈ 20–25) ──────────────────────────────
snap_indices = list(range(20, 26))   # f00020 … f00025
snap_files   = [os.path.join(RUN_DIR, f'field0.f{i:05d}') for i in snap_indices]
snap_files   = [f for f in snap_files if os.path.exists(f)]
print(f'Averaging {len(snap_files)} snapshots: {[os.path.basename(f) for f in snap_files]}')

# ── accumulate mean u over snapshots ────────────────────────────────────────
u_sum = np.zeros_like(y_all)
for fpath in snap_files:
    fld_data = preadnek(fpath, comm)
    # elem.vel shape (3, lz, ly, lx); vel[0] = streamwise u
    u_snap = np.concatenate([e.vel[0].flatten() for e in fld_data.elem])
    u_sum += u_snap
    print(f'  read {os.path.basename(fpath)}, u_max={u_snap.max():.4f}')

u_mean_flat = u_sum / max(1, len(snap_files))

# ── bin by y into ~150 bins spanning [-1, 1] ────────────────────────────────
n_bins = 150
y_edges  = np.linspace(-1.0, 1.0, n_bins + 1)
y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
u_binned  = np.zeros(n_bins)
counts    = np.zeros(n_bins, dtype=int)
for i in range(n_bins):
    mask = (y_all >= y_edges[i]) & (y_all < y_edges[i+1])
    if mask.sum() > 0:
        u_binned[i] = u_mean_flat[mask].mean()
        counts[i]   = mask.sum()

valid = counts > 0
y_prof = y_centers[valid]
u_prof = u_binned[valid]

# ── convert to wall units ────────────────────────────────────────────────────
# Use both walls by folding on |y|; y_wall_dist = 1 - |y|
y_wall = 1.0 - np.abs(y_prof)
yp     = y_wall * Re_tau          # y+
up     = u_prof / u_tau           # u+

# Sort by y+
order = np.argsort(yp)
yp    = yp[order]
up    = up[order]

# ── reference: Reichardt profile ────────────────────────────────────────────
yp_ref = np.logspace(np.log10(0.5), np.log10(Re_tau), 300)
up_ref = reichardt(yp_ref)

# ── reference: viscous sublayer (u+ = y+) ───────────────────────────────────
yp_vis = np.linspace(0, 5, 30)

# ── plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))

ax.semilogx(yp_ref, up_ref, 'k--', lw=1.5, label='Reichardt (model)')
ax.semilogx(yp_vis, yp_vis, 'gray', lw=1.0, ls=':', label='$u^+ = y^+$ (viscous sublayer)')
ax.semilogx(yp, up, 'o', color='#1f4e79', ms=3.5, alpha=0.7,
            label='simulation (time-avg, t=20–25)')

ax.set(xlabel='$y^+$', ylabel='$u^+$',
       title=f'Mean velocity profile  ($Re_\\tau={Re_tau:.0f}$)',
       xlim=(0.5, Re_tau + 20))
ax.legend(fontsize=9)
ax.grid(True, which='both', alpha=0.3)

fig.tight_layout()
out2 = os.path.join(OUT_DIR, 'meanprofile_single_phase.png')
plt.savefig(out2, dpi=150)
plt.close()
print(f'Saved: {out2}')

print('\n--- Summary ---')
print(f'  Snapshots used:  t = 20–25  ({len(snap_files)} files)')
print(f'  u_max (ekin):    {u_max[-5:].mean():.4f} ± {u_max[-5:].std():.4f}  (last 5 rows)')
print(f'  Reichardt U_cl:  {u_cl:.4f}  (expected turbulent centerline)')
print(f'  Laminar U_cl:    1.5000  (Poiseuille)')
print(f'  Conclusion:      {"TURBULENT" if u_max[-5:].mean() < 1.45 else "NOT YET TURBULENT"}')
