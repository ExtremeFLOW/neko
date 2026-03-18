"""
Diagnostics for the two-phase turbulent channel flow case.
Reads ekin.csv from the run directory and produces a diagnostics figure.

Usage:
    python3 analyze_two_phase_channel.py
    python3 analyze_two_phase_channel.py --run channel_test_v4
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--run', default='channel_test_v4',
                    help='Run directory name under /lscratch/sieburgh/simulations/')
args = parser.parse_args()

RUN_DIR = f'/lscratch/sieburgh/simulations/{args.run}'
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# -------------------------------------------------------------------
# Section 1 — ekin.csv time series
# -------------------------------------------------------------------
ekin_path = os.path.join(RUN_DIR, 'ekin.csv')
df = pd.read_csv(
    ekin_path, comment='#',
    names=['t', 'Ekin', 'enst', 'u_max', 'kappa_max', 'kappa_min',
           'kappa_rms', 'Fst_max', 'phi_min', 'phi_max', '_'],
)

print(f"Loaded {len(df)} rows from {ekin_path}")
print(f"  t range:       [{df.t.min():.4f}, {df.t.max():.4f}]")
print(f"  phi_max range: [{df.phi_max.min():.4f}, {df.phi_max.max():.4f}]")
print(f"  phi_min range: [{df.phi_min.min():.4f}, {df.phi_min.max():.4f}]")
print(f"  u_max range:   [{df.u_max.min():.4f}, {df.u_max.max():.4f}]")
print(f"  kappa_rms:     [{df.kappa_rms.min():.4f}, {df.kappa_rms.max():.4f}]")

fig, axes = plt.subplots(2, 2, figsize=(10, 6))

axes[0, 0].plot(df.t, df.u_max)
axes[0, 0].set(ylabel='u_max', xlabel='t', title='Max velocity')

axes[0, 1].plot(df.t, df.Ekin)
axes[0, 1].set(ylabel='E_kin', xlabel='t', title='Kinetic energy')

axes[1, 0].plot(df.t, df.phi_max, label='phi_max')
axes[1, 0].plot(df.t, df.phi_min, label='phi_min')
axes[1, 0].axhline(0.0, color='k', lw=0.5, ls='--')
axes[1, 0].axhline(1.0, color='k', lw=0.5, ls='--')
axes[1, 0].set(ylabel='phi bounds', xlabel='t', title='Phase field bounds')
axes[1, 0].legend()

axes[1, 1].plot(df.t, df.kappa_rms, label='kappa_rms')
axes[1, 1].plot(df.t, df.kappa_max, label='kappa_max', alpha=0.6)
axes[1, 1].set(ylabel='curvature', xlabel='t', title='Curvature')
axes[1, 1].legend()

run_name = os.path.basename(RUN_DIR)
fig.suptitle(f'{run_name}: ekin.csv diagnostics')
plt.tight_layout()

out_path = os.path.join(OUT_DIR, f'diagnostics_ekin_{run_name}.png')
plt.savefig(out_path, dpi=150)
print(f"Saved: {out_path}")
plt.close()

# -------------------------------------------------------------------
# Section 2 — Annotated diagnostics summary
# -------------------------------------------------------------------
print("\n--- Summary ---")
print(f"  kappa_rms(t=0):    {df.kappa_rms.iloc[0]:.2f}  (theoretical 2/R = {2/0.3:.2f})")
print(f"  kappa_rms(t_end):  {df.kappa_rms.iloc[-1]:.2f}")
print(f"  phi_max drop:      {df.phi_max.iloc[0]:.4f} → {df.phi_max.iloc[-1]:.4f}")
print(f"  u_max drop:        {df.u_max.iloc[0]:.4f} → {df.u_max.iloc[-1]:.4f}")
print(f"  Fst_max range:     [{df.Fst_max.min():.4f}, {df.Fst_max.max():.4f}]")

# -------------------------------------------------------------------
# Section 3 — Note on field file visualization (Tier 2, deferred)
# -------------------------------------------------------------------
print("""
Tier 2 (field snapshots) — DEFERRED
=====================================
Field files (field0.f00000 etc.) require pysemtools or Paraview.
pysemtools is not currently installed on egidius.
Options:
  - pip install pysemtools (or from source)
  - Open field0.nek5000 in Paraview if installed
  - Transfer to a machine with pysemtools and use the oscillating_droplet notebook pattern
""")
