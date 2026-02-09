#!/usr/bin/env python3
"""
Simple script to check the initial scalar field values.
"""

import numpy as np
from pathlib import Path
import sys
from mpi4py import MPI

# Import pysemtools
try:
    from pysemtools.io.ppymech.neksuite import preadnek
    from pysemtools.datatypes.field import Field as field_c
    from pysemtools.datatypes.msh import Mesh as msh_c
except ImportError:
    print("Error: pysemtools not found. Please install pysemtools.")
    sys.exit(1)

comm = MPI.COMM_WORLD


def initial_condition_expected(x, y, z, eps):
    """
    Expected initial condition for advecting drop.

    Drop with diameter D=1.0 (radius=0.5) centered at (1.0, 1.0) in 2x2 domain
    Convention: phi = 0 inside droplet, phi = 1 outside (matches Fortran code)
    2D problem: radius only depends on x and y, NOT z
    """
    center_x = 1.0
    center_y = 1.0
    radius = 0.5

    # 2D distance (no z-component!)
    rad = np.sqrt((x - center_x) ** 2 + (y - center_y) ** 2)
    field = 0.5 * (1 + np.tanh((rad - radius) / (2 * eps)))

    return field


def main():
    # Configuration - parse command line arguments if provided
    import sys
    if len(sys.argv) > 1:
        field_file = Path(sys.argv[1])
        epsilon = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    else:
        field_file = Path("../visualization_output/field0.f00000")
        epsilon = 0.05  # From case file

    if not field_file.exists():
        print(f"Error: {field_file} not found")
        print("Make sure you've run the simulation first")
        sys.exit(1)

    print("=" * 80)
    print("Initial Condition Check")
    print("=" * 80)
    print(f"Reading: {field_file}")
    print(f"Expected epsilon: {epsilon}")
    print("=" * 80)

    # Load the field
    data = preadnek(str(field_file), comm)
    msh = msh_c(comm, data=data)
    fld = field_c(comm, data=data)

    # Get the scalar field
    phi_actual = fld.fields["scal"][0]

    # Compute expected initial condition
    phi_expected = initial_condition_expected(msh.x, msh.y, msh.z, epsilon)

    # Print statistics for actual field
    print("\nActual scalar field (from simulation):")
    print(f"  Min:    {np.min(phi_actual):.6f}")
    print(f"  Max:    {np.max(phi_actual):.6f}")
    print(f"  Mean:   {np.mean(phi_actual):.6f}")
    print(f"  Median: {np.median(phi_actual):.6f}")

    # Count points in different ranges
    n_below_0 = np.sum(phi_actual < 0)
    n_above_1 = np.sum(phi_actual > 1)
    n_in_range = np.sum((phi_actual >= 0) & (phi_actual <= 1))
    n_total = len(phi_actual)

    print(f"\nBound check:")
    print(f"  Points < 0:   {n_below_0} ({100*n_below_0/n_total:.2f}%)")
    print(f"  Points [0,1]: {n_in_range} ({100*n_in_range/n_total:.2f}%)")
    print(f"  Points > 1:   {n_above_1} ({100*n_above_1/n_total:.2f}%)")

    # Print statistics for expected field
    print("\n" + "-" * 80)
    print("Expected scalar field (analytical):")
    print(f"  Min:    {np.min(phi_expected):.6f}")
    print(f"  Max:    {np.max(phi_expected):.6f}")
    print(f"  Mean:   {np.mean(phi_expected):.6f}")
    print(f"  Median: {np.median(phi_expected):.6f}")

    # Compare actual vs expected
    print("\n" + "-" * 80)
    print("Comparison (actual vs expected):")
    diff = phi_actual - phi_expected
    print(f"  Max difference: {np.max(np.abs(diff)):.6e}")
    print(f"  RMS difference: {np.sqrt(np.mean(diff**2)):.6e}")

    if np.max(np.abs(diff)) < 1e-10:
        print("\n✓ Initial condition matches analytical formula!")
    else:
        print("\n⚠ Initial condition differs from analytical formula")
        print("  This might indicate an issue with the Fortran initial_conditions subroutine")

    # Check domain extents
    print("\n" + "-" * 80)
    print("Domain extents:")
    print(f"  x: [{np.min(msh.x):.3f}, {np.max(msh.x):.3f}]")
    print(f"  y: [{np.min(msh.y):.3f}, {np.max(msh.y):.3f}]")
    print(f"  z: [{np.min(msh.z):.3f}, {np.max(msh.z):.3f}]")

    # Check droplet center location (where phi ≈ 0.5)
    interface_mask = np.abs(phi_actual - 0.5) < 0.1
    if np.sum(interface_mask) > 0:
        center_x_approx = np.mean(msh.x[interface_mask])
        center_y_approx = np.mean(msh.y[interface_mask])
        print(f"\nApproximate droplet center (where φ ≈ 0.5):")
        print(f"  Center: ({center_x_approx:.3f}, {center_y_approx:.3f})")
        print(f"  Expected: (1.0, 1.0)")

    print("\n" + "=" * 80)
    print("Check complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
