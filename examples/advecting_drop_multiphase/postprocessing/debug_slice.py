#!/usr/bin/env python3
"""Debug script to check actual values in 2D slice."""

import numpy as np
from pathlib import Path
from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.field import Field as field_c
from pysemtools.datatypes.msh import Mesh as msh_c

comm = MPI.COMM_WORLD

# Load initial field
field_file = Path("../visualization_output/field0.f00000")
data = preadnek(str(field_file), comm)
msh = msh_c(comm, data=data)
fld = field_c(comm, data=data)

phi = fld.fields["scal"][0]

# Extract 2D slice at z=0 (same as animation)
z_tol = 0.1
z_mask = np.abs(msh.z) < z_tol

x_2d = msh.x[z_mask]
y_2d = msh.y[z_mask]
phi_2d = phi[z_mask]

print("=" * 80)
print("2D Slice Analysis (|z| < 0.1)")
print("=" * 80)
print(f"Number of points in slice: {len(phi_2d)}")
print(f"Min φ in slice:  {np.min(phi_2d):.6f}")
print(f"Max φ in slice:  {np.max(phi_2d):.6f}")
print(f"Mean φ in slice: {np.mean(phi_2d):.6f}")
print()

# Find the point closest to droplet center (1, 1, 0)
center_x, center_y = 1.0, 1.0
distances = np.sqrt((x_2d - center_x)**2 + (y_2d - center_y)**2)
closest_idx = np.argmin(distances)

z_2d = msh.z[z_mask]
print("Point closest to droplet center (1, 1):")
print(f"  Location: ({x_2d[closest_idx]:.6f}, {y_2d[closest_idx]:.6f}, {z_2d[closest_idx]:.6f})")
print(f"  φ value:  {phi_2d[closest_idx]:.6f}")
print(f"  Distance from center: {distances[closest_idx]:.6f}")
print()

# Check what z values are in the slice
print(f"Z-coordinates in slice range: [{np.min(z_2d):.6f}, {np.max(z_2d):.6f}]")
print(f"Number of points with z ≈ 0 (|z| < 0.01): {np.sum(np.abs(z_2d) < 0.01)}")
print()

# Check points at different radii
print("φ values at different radii from center (1,1):")
for radius in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
    mask = (distances >= radius - 0.05) & (distances <= radius + 0.05)
    if np.sum(mask) > 0:
        phi_avg = np.mean(phi_2d[mask])
        print(f"  r = {radius:.1f}: φ = {phi_avg:.6f} ({np.sum(mask)} points)")

print()
print("=" * 80)

# Check full 3D field for comparison
print("Full 3D Field (for comparison):")
print(f"Min φ:  {np.min(phi):.6f}")
print(f"Max φ:  {np.max(phi):.6f}")
print("=" * 80)
