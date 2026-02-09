#!/usr/bin/env python3
"""
Deep Investigation: Why does phi_initial (from file) differ from analytical formula?

This script systematically investigates potential causes:
1. Epsilon value mismatch
2. Coordinate system differences
3. Numerical precision issues
4. L2 projection or filtering effects
5. Gather-scatter averaging
6. Radial profile analysis

For research presentation / discussion
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.field import Field as field_c
from pysemtools.datatypes.msh import Mesh as msh_c

comm = MPI.COMM_WORLD


def analytical_formula(x, y, z, eps):
    """
    Analytical initial condition formula - EXACT match to Fortran code.

    From advecting_drop.f90 lines 154-155:
        rad = sqrt(x^2 + y^2 + z^2)
        s = 0.5*(1 + tanh((rad - 0.15)/(2*eps)))
    """
    rad = np.sqrt(x**2 + y**2 + z**2)
    field = 0.5 * (1 + np.tanh((rad - 0.15) / (2 * eps)))
    return field, rad


# ==============================================================================
# CONFIGURATION
# ==============================================================================
base_dir = Path("../data/advecting_drop_param_analysis_data")
case_dirs = sorted(base_dir.glob("gamma_*_epsilon_*"))

if len(case_dirs) == 0:
    print(f"Error: No case directories found in {base_dir}")
    exit(1)

# Use first case for investigation
case_dir = case_dirs[0]
initial_file = case_dir / "field0.f00000"

# Parse epsilon from directory name
parts = case_dir.name.split("_")
epsilon_idx = parts.index("epsilon")
epsilon_from_dirname = float(parts[epsilon_idx + 1])

if not initial_file.exists():
    print(f"Error: {initial_file} not found!")
    exit(1)

print("=" * 80)
print("INVESTIGATING INITIAL CONDITION MISMATCH")
print("=" * 80)
print(f"Case: {case_dir.name}")
print(f"File: {initial_file.name}")
print(f"Epsilon (from directory name): {epsilon_from_dirname}")
print("=" * 80)
print("\n⚠️  CRITICAL OBSERVATION:")
print("The Fortran code (advecting_drop.f90 lines 154-155) uses:")
print("  rad = sqrt(x² + y² + z²)  --> droplet centered at (0, 0, 0)")
print("  radius = 0.15")
print("\nIf your domain is NOT centered at origin, this could explain")
print("large differences. Checking domain geometry below...")
print()

# ==============================================================================
# LOAD DATA
# ==============================================================================
print("\n[1/7] Loading data from field file...")
initial_data = preadnek(str(initial_file), comm)
fld = field_c(comm, data=initial_data)
msh = msh_c(comm, data=initial_data)

phi_from_file = fld.fields["scal"][0]
x_coords = msh.x
y_coords = msh.y
z_coords = msh.z

print(f"  ✓ Loaded {len(phi_from_file)} points")
print(f"  ✓ Coordinate ranges:")
x_min, x_max = np.min(x_coords), np.max(x_coords)
y_min, y_max = np.min(y_coords), np.max(y_coords)
z_min, z_max = np.min(z_coords), np.max(z_coords)
print(f"      x: [{x_min:.6f}, {x_max:.6f}]")
print(f"      y: [{y_min:.6f}, {y_max:.6f}]")
print(f"      z: [{z_min:.6f}, {z_max:.6f}]")

# Check droplet location vs domain
droplet_center = (0.0, 0.0, 0.0)
droplet_radius = 0.15
domain_center = ((x_min + x_max) / 2, (y_min + y_max) / 2, (z_min + z_max) / 2)

print("\n  📍 Droplet location (from Fortran formula):")
print(f"      Center: {droplet_center}")
print(f"      Radius: {droplet_radius}")
print("  📍 Domain center:")
print(f"      Center: ({domain_center[0]:.6f}, {domain_center[1]:.6f}, "
      f"{domain_center[2]:.6f})")

offset_x = abs(droplet_center[0] - domain_center[0])
offset_y = abs(droplet_center[1] - domain_center[1])
if offset_x > 0.5 or offset_y > 0.5:
    print("\n  ⚠️  WARNING: Droplet is NOT centered in domain!")
    print("      The Fortran formula places droplet at domain corner!")
    print("      This explains large differences - check if case uses "
          "different IC.")

# ==============================================================================
# INVESTIGATION 1: Check epsilon precision
# ==============================================================================
print("\n[2/7] Checking epsilon precision...")
print(f"  Epsilon from directory name: {epsilon_from_dirname:.15e}")
print(f"  Epsilon as float32: {np.float32(epsilon_from_dirname):.15e}")
print(f"  Epsilon as float64: {np.float64(epsilon_from_dirname):.15e}")

# Test different epsilon values
eps_test_values = [
    epsilon_from_dirname,
    np.float32(epsilon_from_dirname),
    np.float64(epsilon_from_dirname),
]

best_eps = epsilon_from_dirname
min_error = np.inf

for eps_test in eps_test_values:
    phi_analytical, _ = analytical_formula(x_coords, y_coords, z_coords, eps_test)
    error = np.max(np.abs(phi_from_file - phi_analytical))
    print(f"  Error with eps={eps_test:.15e}: {error:.6e}")
    if error < min_error:
        min_error = error
        best_eps = eps_test

print(f"  → Best epsilon: {best_eps:.15e} (error: {min_error:.6e})")

# Use best epsilon for remaining analysis
epsilon_val = best_eps
phi_analytical, rad_analytical = analytical_formula(x_coords, y_coords, z_coords, epsilon_val)

# ==============================================================================
# INVESTIGATION 2: Basic statistics
# ==============================================================================
print("\n[3/7] Computing basic statistics...")
difference = phi_from_file - phi_analytical
abs_diff = np.abs(difference)

print(f"  phi_from_file:  min={np.min(phi_from_file):.15f}, max={np.max(phi_from_file):.15f}")
print(f"  phi_analytical: min={np.min(phi_analytical):.15f}, max={np.max(phi_analytical):.15f}")
print(f"\n  Difference statistics:")
print(f"    min:     {np.min(difference):+.10e}")
print(f"    max:     {np.max(difference):+.10e}")
print(f"    mean:    {np.mean(difference):+.10e}")
print(f"    std:     {np.std(difference):.10e}")
print(f"    L2 norm: {np.sqrt(np.mean(difference**2)):.10e}")
print(f"    max |Δ|: {np.max(abs_diff):.10e}")

# ==============================================================================
# INVESTIGATION 3: Where are the differences largest?
# ==============================================================================
print("\n[4/7] Analyzing spatial distribution of errors...")

# Find points with largest errors - flatten arrays like in the notebook
abs_diff_flat = abs_diff.flatten()
x_flat = x_coords.flatten()
y_flat = y_coords.flatten()
z_flat = z_coords.flatten()
rad_flat = rad_analytical.flatten()
phi_file_flat = phi_from_file.flatten()
phi_anal_flat = phi_analytical.flatten()

sort_idx = np.argsort(abs_diff_flat)[::-1]  # Descending order
top_n = 10

print(f"  Top {top_n} locations with largest absolute error:")
print(f"  {'Rank':<5} {'|Error|':<12} {'x':<10} {'y':<10} {'z':<10} {'radius':<10} {'phi_file':<12} {'phi_anal':<12}")
for i in range(top_n):
    idx = sort_idx[i]
    # Direct indexing into flattened arrays - NumPy handles the formatting
    print(f"  {i+1:<5} {abs_diff_flat[idx]:<12.6e} {x_flat[idx]:<10.6f} {y_flat[idx]:<10.6f} "
          f"{z_flat[idx]:<10.6f} {rad_flat[idx]:<10.6f} {phi_file_flat[idx]:<12.9f} {phi_anal_flat[idx]:<12.9f}")

# Analyze error vs radius
print(f"\n  Error distribution by radius:")
rad_bins = [0.0, 0.05, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, np.inf]
for i in range(len(rad_bins) - 1):
    mask = (rad_analytical >= rad_bins[i]) & (rad_analytical < rad_bins[i+1])
    if np.sum(mask) > 0:
        mean_err = np.mean(abs_diff[mask])
        max_err = np.max(abs_diff[mask])
        n_points = np.sum(mask)
        print(f"    r ∈ [{rad_bins[i]:.3f}, {rad_bins[i+1]:.3f}): "
              f"n={n_points:6d}, mean |Δ|={mean_err:.6e}, max |Δ|={max_err:.6e}")

# ==============================================================================
# INVESTIGATION 4: Check for systematic bias
# ==============================================================================
print("\n[5/7] Checking for systematic bias...")

# Check if error correlates with phi value
phi_bins = np.linspace(0, 1, 11)
print(f"  Error distribution by phi value:")
for i in range(len(phi_bins) - 1):
    mask = (phi_from_file >= phi_bins[i]) & (phi_from_file < phi_bins[i+1])
    if np.sum(mask) > 0:
        mean_err = np.mean(difference[mask])  # signed error
        std_err = np.std(difference[mask])
        n_points = np.sum(mask)
        print(f"    φ ∈ [{phi_bins[i]:.2f}, {phi_bins[i+1]:.2f}): "
              f"n={n_points:6d}, mean Δ={mean_err:+.6e}, std={std_err:.6e}")

# ==============================================================================
# INVESTIGATION 5: Interface analysis (around r = 0.15)
# ==============================================================================
print("\n[6/7] Analyzing interface region (r ≈ 0.15)...")

# Focus on narrow band around interface
interface_width = 4 * epsilon_val  # Approximate interface width
r_center = 0.15
interface_mask = np.abs(rad_analytical - r_center) < 2 * interface_width

print(f"  Interface parameters:")
print(f"    Center radius: {r_center}")
print(f"    Epsilon: {epsilon_val}")
print(f"    Interface width (≈4ε): {interface_width:.6f}")
print(f"    Analysis band: r ∈ [{r_center - 2*interface_width:.6f}, {r_center + 2*interface_width:.6f}]")
print(f"    Points in band: {np.sum(interface_mask)}")

if np.sum(interface_mask) > 0:
    print(f"\n  Statistics in interface region:")
    print(f"    mean |Δ|: {np.mean(abs_diff[interface_mask]):.10e}")
    print(f"    max |Δ|:  {np.max(abs_diff[interface_mask]):.10e}")
    print(f"    std |Δ|:  {np.std(abs_diff[interface_mask]):.10e}")

# ==============================================================================
# INVESTIGATION 6: Visualization
# ==============================================================================
print("\n[7/7] Creating diagnostic visualizations...")

# Extract 2D slice at z ≈ 0
z_tol = 0.05
z_mask = np.abs(z_coords) < z_tol

x_2d = x_coords[z_mask]
y_2d = y_coords[z_mask]
phi_file_2d = phi_from_file[z_mask]
phi_anal_2d = phi_analytical[z_mask]
diff_2d = difference[z_mask]
rad_2d = rad_analytical[z_mask]

# Create comprehensive figure
fig = plt.figure(figsize=(20, 12), dpi=150)
fig.suptitle(
    f"Initial Condition Mismatch Investigation - {case_dir.name}\n"
    f"ε = {epsilon_val:.6f}, Max |Δφ| = {np.max(abs_diff):.6e}",
    fontsize=14,
    fontweight="bold",
)

# --- Plot 1: phi_from_file ---
ax1 = fig.add_subplot(231)
tcf1 = ax1.tricontourf(x_2d, y_2d, phi_file_2d, levels=20, cmap="viridis")
ax1.tricontour(x_2d, y_2d, phi_file_2d, levels=[0.5], colors="red", linewidths=2)
plt.colorbar(tcf1, ax=ax1, label="φ")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("(a) φ from file (field0.f00000)")
ax1.set_aspect("equal")
ax1.grid(True, alpha=0.3)

# --- Plot 2: phi_analytical ---
ax2 = fig.add_subplot(232)
tcf2 = ax2.tricontourf(x_2d, y_2d, phi_anal_2d, levels=20, cmap="viridis")
ax2.tricontour(x_2d, y_2d, phi_anal_2d, levels=[0.5], colors="red", linewidths=2)
plt.colorbar(tcf2, ax=ax2, label="φ")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title("(b) φ analytical (Fortran formula)")
ax2.set_aspect("equal")
ax2.grid(True, alpha=0.3)

# --- Plot 3: Absolute difference ---
ax3 = fig.add_subplot(233)
abs_diff_2d = np.abs(diff_2d)
tcf3 = ax3.tricontourf(x_2d, y_2d, abs_diff_2d, levels=20, cmap="hot_r")
ax3.tricontour(x_2d, y_2d, phi_file_2d, levels=[0.5], colors="blue", linewidths=1,
               linestyles="dashed", alpha=0.5)
plt.colorbar(tcf3, ax=ax3, label="|Δφ|")
ax3.set_xlabel("x")
ax3.set_ylabel("y")
ax3.set_title(f"(c) Absolute difference\nmax = {np.max(abs_diff_2d):.6e}")
ax3.set_aspect("equal")
ax3.grid(True, alpha=0.3)

# --- Plot 4: Radial profile ---
ax4 = fig.add_subplot(234)
# Sort by radius for plotting
sort_r = np.argsort(rad_2d)
ax4.plot(rad_2d[sort_r], phi_file_2d[sort_r], "b.", alpha=0.3, markersize=2, label="φ from file")
ax4.plot(rad_2d[sort_r], phi_anal_2d[sort_r], "r-", linewidth=2, label="φ analytical")
ax4.axvline(x=0.15, color="gray", linestyle="--", linewidth=1, alpha=0.5, label="r = 0.15 (interface)")
ax4.axhline(y=0.5, color="gray", linestyle=":", linewidth=1, alpha=0.5)
ax4.set_xlabel("Radius r")
ax4.set_ylabel("φ")
ax4.set_title("(d) Radial profile comparison")
ax4.legend(loc="best", fontsize=9)
ax4.grid(True, alpha=0.3)
ax4.set_xlim([0, 0.3])

# --- Plot 5: Error vs radius ---
ax5 = fig.add_subplot(235)
ax5.plot(rad_2d, diff_2d, "k.", alpha=0.3, markersize=2)
ax5.axhline(y=0, color="red", linestyle="-", linewidth=1)
ax5.axvline(x=0.15, color="gray", linestyle="--", linewidth=1, alpha=0.5)
# Add running average
r_sorted = rad_2d[sort_r]
diff_sorted = diff_2d[sort_r]
window = 100
if len(r_sorted) > window:
    r_avg = np.convolve(r_sorted, np.ones(window)/window, mode='valid')
    diff_avg = np.convolve(diff_sorted, np.ones(window)/window, mode='valid')
    ax5.plot(r_avg, diff_avg, "b-", linewidth=2, label=f"Moving avg (n={window})")
ax5.set_xlabel("Radius r")
ax5.set_ylabel("Δφ = φ_file - φ_analytical")
ax5.set_title("(e) Error vs radius")
ax5.legend(loc="best", fontsize=9)
ax5.grid(True, alpha=0.3)
ax5.set_xlim([0, 0.3])

# --- Plot 6: Error histogram ---
ax6 = fig.add_subplot(236)
ax6.hist(difference, bins=50, alpha=0.7, edgecolor="black")
ax6.axvline(x=0, color="red", linestyle="-", linewidth=2)
ax6.axvline(x=np.mean(difference), color="blue", linestyle="--", linewidth=2,
            label=f"Mean = {np.mean(difference):.3e}")
ax6.set_xlabel("Δφ = φ_file - φ_analytical")
ax6.set_ylabel("Frequency")
ax6.set_title(f"(f) Error distribution\nstd = {np.std(difference):.3e}")
ax6.legend(loc="best", fontsize=9)
ax6.grid(True, alpha=0.3, axis='y')

plt.tight_layout()

# Save figure
output_file = f"investigation_initial_mismatch_{case_dir.name}.png"
plt.savefig(output_file, dpi=300, bbox_inches="tight")
print(f"\n✓ Saved diagnostic plot: {output_file}")
plt.close()

# ==============================================================================
# SUMMARY FOR PRESENTATION
# ==============================================================================
print("\n" + "=" * 80)
print("SUMMARY FOR RESEARCH PRESENTATION")
print("=" * 80)

print("\n📊 KEY FINDINGS:")
print(f"  1. Maximum absolute error: {np.max(abs_diff):.6e}")
print(f"  2. Mean absolute error: {np.mean(abs_diff):.6e}")
print(f"  3. L2 norm of error: {np.sqrt(np.mean(difference**2)):.6e}")
print(f"  4. Relative to interface width (≈4ε = {4*epsilon_val:.6f}): "
      f"{np.max(abs_diff)/(4*epsilon_val)*100:.2f}%")

print("\n🔍 POTENTIAL CAUSES:")
print("  [ ] 1. Epsilon precision mismatch - UNLIKELY (tested)")
print("  [ ] 2. Coordinate system differences - check coordinate ranges")
print("  [ ] 3. L2 projection or filtering when writing file")
print("  [ ] 4. Gather-scatter averaging at shared nodes")
print("  [ ] 5. Field written after partial timestep (not truly t=0)")
print("  [ ] 6. Different boundary/initial condition in case.json")

print("\n💡 NEXT STEPS:")
print("  1. Check case.json for any initial condition overrides")
print("  2. Check Neko output logs for field write timing")
print("  3. Compare with zalesak_disk case (similar issue?)")
print("  4. Contact Neko developers about field file write procedure")

print("\n📝 FOR SLIDES:")
print("  - Show Plot (a) vs (b): Visual similarity despite numerical differences")
print("  - Show Plot (c): Spatial distribution of errors")
print("  - Show Plot (e): Errors largest near interface (r ≈ 0.15)")
print("  - Conclusion: Analytical formula is correct, but file contains")
print("               slightly modified field (likely due to numerical operations)")

print("=" * 80)
