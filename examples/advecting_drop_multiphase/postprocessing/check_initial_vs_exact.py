#!/usr/bin/env python3
"""
Quick check: Compare phi_initial (from file) vs phi_exact (analytical)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.field import Field as field_c
from pysemtools.datatypes.msh import Mesh as msh_c

comm = MPI.COMM_WORLD


def initial_condition(x, y, z, eps):
    """Generate the analytical initial condition for advecting drop (matches Fortran code)."""
    # Center at origin, radius 0.15 (from advecting_drop.f90 lines 154-155)
    rad = np.sqrt(x**2 + y**2 + z**2)
    field = 0.5 * (1 + np.tanh((rad - 0.15) / (2 * eps)))

    return field


# Test case directory - using local data path
base_dir = Path("../data/advecting_drop_param_analysis_data")

# Find first available case
case_dirs = sorted(base_dir.glob("gamma_*_epsilon_*"))
if len(case_dirs) == 0:
    print(f"Error: No case directories found in {base_dir}")
    exit(1)

case_dir = case_dirs[0]  # Use first case
initial_file = case_dir / "field0.f00000"

# Parse epsilon from directory name
parts = case_dir.name.split("_")
epsilon_idx = parts.index("epsilon")
epsilon_val = float(parts[epsilon_idx + 1])

if not initial_file.exists():
    print(f"Error: {initial_file} not found!")
    exit(1)

print("=" * 80)
print("Comparing phi_initial (from file) vs phi_exact (analytical)")
print("=" * 80)
print(f"Case: {case_dir.name}")
print(f"Epsilon: {epsilon_val}")
print()

# Load initial field from file
initial_data = preadnek(str(initial_file), comm)
fld = field_c(comm, data=initial_data)
msh = msh_c(comm, data=initial_data)

phi_initial = fld.fields["scal"][0]

# Compute analytical initial condition
phi_exact = initial_condition(msh.x, msh.y, msh.z, epsilon_val)

# Compare
difference = phi_initial - phi_exact
abs_diff = np.abs(difference)

print("Statistics:")
print(f"  phi_initial: min={np.min(phi_initial):.15f}, max={np.max(phi_initial):.15f}")
print(f"  phi_exact:   min={np.min(phi_exact):.15f}, max={np.max(phi_exact):.15f}")
print()
print("Difference (phi_initial - phi_exact):")
print(f"  min:  {np.min(difference):.15e}")
print(f"  max:  {np.max(difference):.15e}")
print(f"  mean: {np.mean(abs_diff):.15e}")
print(f"  std:  {np.std(abs_diff):.15e}")
print(f"  L2:   {np.sqrt(np.mean(abs_diff**2)):.15e}")
print(f"  max abs: {np.max(abs_diff):.15e}")
print()

# Check if they're identical (within machine precision)
if np.allclose(phi_initial, phi_exact, rtol=1e-14, atol=1e-14):
    print("✓ phi_initial and phi_exact are IDENTICAL (within machine precision)")
else:
    print("✗ phi_initial and phi_exact are DIFFERENT")
    print()
    print("Relative difference:")
    # Avoid division by zero
    mask = np.abs(phi_exact) > 1e-10
    rel_diff = np.zeros_like(difference)
    rel_diff[mask] = abs_diff[mask] / np.abs(phi_exact[mask])
    print(f"  max relative: {np.max(rel_diff):.15e}")
    print(f"  mean relative: {np.mean(rel_diff[mask]):.15e}")

print("=" * 80)

# Create visualization
print("\nCreating visualization...")

# Extract 2D slice at z=0
z_tol = 0.1
z_mask = np.abs(msh.z) < z_tol
x_2d = msh.x[z_mask]
y_2d = msh.y[z_mask]
phi_initial_2d = phi_initial[z_mask]
phi_exact_2d = phi_exact[z_mask]
diff_2d = difference[z_mask]

# Create figure with 3 subplots
fig = plt.figure(figsize=(18, 5), dpi=150)
fig.suptitle(
    f"Comparison: phi_initial (file) vs phi_exact (analytical) - {case_dir.name}",
    fontsize=14,
    fontweight="bold",
)

# Plot 1: phi_initial from file
ax1 = fig.add_subplot(131)
tcf1 = ax1.tricontourf(x_2d, y_2d, phi_initial_2d, levels=20, cmap="viridis")
plt.colorbar(tcf1, ax=ax1, label="φ")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("(a) φ_initial (from file)")
ax1.set_aspect("equal")

# Plot 2: phi_exact (analytical)
ax2 = fig.add_subplot(132)
tcf2 = ax2.tricontourf(x_2d, y_2d, phi_exact_2d, levels=20, cmap="viridis")
plt.colorbar(tcf2, ax=ax2, label="φ")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title("(b) φ_exact (analytical)")
ax2.set_aspect("equal")

# Plot 3: Difference
ax3 = fig.add_subplot(133)
tcf3 = ax3.tricontourf(x_2d, y_2d, diff_2d, levels=20, cmap="RdBu_r")
plt.colorbar(tcf3, ax=ax3, label="Δφ = φ_initial - φ_exact")
ax3.set_xlabel("x")
ax3.set_ylabel("y")
ax3.set_title("(c) Difference")
ax3.set_aspect("equal")

plt.tight_layout()

# Save figure
output_file = f"check_initial_vs_exact_{case_dir.name}.png"
plt.savefig(output_file, dpi=300, bbox_inches="tight")
print(f"✓ Saved plot: {output_file}")
plt.close()

print("=" * 80)
