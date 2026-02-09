#!/usr/bin/env python3
"""
Reference-Style Analysis for Advecting Drop Parameter Sweep

Reproduces the analysis from the reference paper:
1. Contour plots of φ=0.5 comparing initial vs final
2. 1D profiles of φ along a line (e.g., y=1.0)
3. Visual comparison of initial and final solutions
"""

import numpy as np
import matplotlib.pyplot as plt
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


def initial_condition(x, y, z, eps):
    """Generate the initial condition for advecting drop."""
    center_x = 1.0
    center_y = 1.0
    center_z = 0.0
    radius = 0.5

    rad = np.sqrt((x - center_x) ** 2 + (y - center_y) ** 2 + (z - center_z) ** 2)
    field = 0.5 * (1 + np.tanh((rad - radius) / (2 * eps)))

    return field


def plot_case_comparison(case_dir, output_dir, gamma_val, epsilon_val, u_max, delta_x):
    """
    Create reference-style plots for a single case.

    Args:
        case_dir: Path to case directory
        output_dir: Path to save plots
        gamma_val: Value of gamma for this case
        epsilon_val: Value of epsilon for this case
        u_max: Maximum velocity for non-dimensionalization
        delta_x: Grid spacing for non-dimensionalization
    """
    case_name = case_dir.name
    print(f"\nAnalyzing {case_name}...")

    # Compute non-dimensional parameters
    gamma_nd = gamma_val / u_max
    epsilon_nd = epsilon_val / delta_x

    # Load initial and final fields
    initial_file = case_dir / "field0.f00000"

    # Find last field file
    field_files = sorted(case_dir.glob("field0.f*"))
    if len(field_files) == 0:
        print(f"  ⚠ Skipping {case_name}: no field files found")
        return
    final_file = field_files[-1]

    if not initial_file.exists() or not final_file.exists():
        print(f"  ⚠ Skipping {case_name}: missing files")
        return

    print(f"  Initial: {initial_file.name}")
    print(f"  Final: {final_file.name}")

    try:
        # Load initial field from file
        initial_data = preadnek(str(initial_file), comm)
        fld_initial = field_c(comm, data=initial_data)
        msh = msh_c(comm, data=initial_data)
        phi_initial = fld_initial.fields["scal"][0]

        # Load final field
        final_data = preadnek(str(final_file), comm)
        fld_final = field_c(comm, data=final_data)
        phi_final = fld_final.fields["scal"][0]

        # Create figure with 3 subplots
        fig = plt.figure(figsize=(18, 5), dpi=150)
        fig.suptitle(
            f"Advecting Drop: γ/u_max = {gamma_nd:.0e}, ε/Δx = {epsilon_nd:.2f}",
            fontsize=14,
            fontweight="bold",
            y=1.02,
        )

        # --- Plot 1: Contour comparison (2D slice at z=0) ---
        ax1 = fig.add_subplot(131)

        # Find indices near z=0
        z_tol = 0.1
        z_mask = np.abs(msh.z) < z_tol

        # Get 2D data
        x_2d = msh.x[z_mask]
        y_2d = msh.y[z_mask]
        phi_initial_2d = phi_initial[z_mask]
        phi_final_2d = phi_final[z_mask]

        # Plot φ=0.5 contours
        ax1.tricontour(
            x_2d,
            y_2d,
            phi_initial_2d,
            levels=[0.5],
            colors="black",
            linewidths=2,
            linestyles="solid",
            label="Initial (t=0)",
        )
        ax1.tricontour(
            x_2d,
            y_2d,
            phi_final_2d,
            levels=[0.5],
            colors="red",
            linewidths=2,
            linestyles="dashed",
            label="Final (t=T)",
        )

        ax1.set_xlabel("x", fontsize=12)
        ax1.set_ylabel("y", fontsize=12)
        ax1.set_title(
            "(a) Droplet Interface Position\nφ=0.5 contour comparison",
            fontsize=11,
            fontweight="bold",
        )
        ax1.set_aspect("equal")
        ax1.legend(loc="lower left", fontsize=9)
        ax1.grid(True, alpha=0.3)

        # --- Plot 2: 1D Profile at y=1.0 ---
        ax2 = fig.add_subplot(132)

        # Extract profile at y ≈ 1.0
        y_tol = 0.05
        y_mask = (np.abs(msh.y - 1.0) < y_tol) & z_mask

        if np.sum(y_mask) > 10:
            x_profile = msh.x[y_mask]
            phi_initial_profile = phi_initial[y_mask]
            phi_final_profile = phi_final[y_mask]

            # Sort by x for plotting
            sort_idx = np.argsort(x_profile)

            ax2.plot(
                x_profile[sort_idx],
                phi_initial_profile[sort_idx],
                "k-",
                linewidth=2,
                label="Initial (t=0)",
            )
            ax2.plot(
                x_profile[sort_idx],
                phi_final_profile[sort_idx],
                "r--",
                linewidth=2,
                label="Final (t=T)",
            )
            ax2.axhline(y=0.5, color="gray", linestyle=":", linewidth=1, alpha=0.5)
            ax2.axhline(y=0.0, color="gray", linestyle="-", linewidth=0.5, alpha=0.3)
            ax2.axhline(y=1.0, color="gray", linestyle="-", linewidth=0.5, alpha=0.3)

            ax2.set_xlabel("x", fontsize=12)
            ax2.set_ylabel("φ (scalar field)", fontsize=12)
            ax2.set_title(
                "(b) 1D Profile Through Droplet\nCut at y=1.0",
                fontsize=11,
                fontweight="bold",
            )
            ax2.legend(loc="best", fontsize=9)
            ax2.grid(True, alpha=0.3)
            ax2.set_ylim([-0.1, 1.1])

            # Add min/max values text
            phi_initial_min = np.min(phi_initial)
            phi_initial_max = np.max(phi_initial)
            phi_final_min = np.min(phi_final)
            phi_final_max = np.max(phi_final)

            stats_text = (
                f"Initial: min={phi_initial_min:.10f}, max={phi_initial_max:.10f}\n"
                f"Final: min={phi_final_min:.10f}, max={phi_final_max:.10f}"
            )

            # Print stats to console instead of overlaying on plot
            print(f"  {stats_text.replace(chr(10), ' | ')}")

        # --- Plot 3: 2D field comparison (filled contours) ---
        ax3 = fig.add_subplot(133)

        # Plot absolute difference field
        diff_2d = np.abs(phi_final_2d - phi_initial_2d)

        # Use fixed levels for consistent colorbar across all cases
        # Range [0, 0.1] since errors > 0.1 are already considered poor
        error_levels = np.linspace(0, 0.1, 21)
        tcf = ax3.tricontourf(
            x_2d, y_2d, diff_2d, levels=error_levels, cmap="hot_r", extend="max"
        )

        plt.colorbar(tcf, ax=ax3, label="|Δφ| = |φ_final - φ_initial|")
        ax3.set_xlabel("x", fontsize=12)
        ax3.set_ylabel("y", fontsize=12)
        ax3.set_title(
            "(c) Error Field\nChange from initial condition",
            fontsize=11,
            fontweight="bold",
        )
        ax3.set_aspect("equal")

        plt.tight_layout()

        # Add figure caption below plots
        fig.text(
            0.5,
            -0.02,
            "Note: Overlapping contours indicate good shape preservation. "
            "Bright regions in error field show numerical diffusion.",
            ha="center",
            fontsize=9,
            style="italic",
        )

        # Save figure
        output_file = output_dir / f"{case_name}_reference_analysis.png"
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"  ✓ Saved: {output_file.name}")
        plt.close()

    except Exception as e:
        print(f"  ✗ Error processing {case_name}: {e}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Reference-style analysis for advecting drop simulation"
    )
    parser.add_argument(
        "--case-dir",
        type=str,
        default=None,
        help="Single case directory to analyze (e.g., ../gamma_0.01_epsilon_0.02)",
    )
    parser.add_argument(
        "--base-dir",
        type=str,
        default="../advecting_drop_param_analysis_data",
        help="Base directory for parameter sweep (default: ../advecting_drop_param_analysis_data)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="reference_style_plots",
        help="Output directory for plots (default: reference_style_plots)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    # Nondimensionalization constants
    u_max = 1.0  # Maximum velocity
    delta_x = 0.00357  # Average GLL spacing for 80x80 mesh on 2x2 domain

    print("=" * 80)
    print("Reference-Style Analysis for Advecting Drop")
    print("=" * 80)
    print(f"Output directory: {output_dir.resolve()}")
    print(f"Non-dimensionalization: u_max = {u_max}, delta_x = {delta_x}")
    print("=" * 80)

    # Determine which cases to analyze
    if args.case_dir:
        # Analyze single case
        case_dir = Path(args.case_dir)
        if not case_dir.exists():
            print(f"Error: Case directory not found: {case_dir}")
            sys.exit(1)

        case_dirs = [case_dir]
        print(f"\nAnalyzing single case: {case_dir.name}")
    else:
        # Analyze parameter sweep
        base_dir = Path(args.base_dir)
        case_dirs = sorted(base_dir.glob("gamma_*_epsilon_*"))

        if len(case_dirs) == 0:
            print(f"Error: No case directories found in {base_dir}")
            sys.exit(1)

        print(f"\nFound {len(case_dirs)} case directories")

    # Process each case
    for case_dir in case_dirs:
        # Parse gamma and epsilon values from directory name
        parts = case_dir.name.split("_")
        try:
            gamma_idx = parts.index("gamma")
            epsilon_idx = parts.index("epsilon")
            gamma_val = float(parts[gamma_idx + 1])
            epsilon_val = float(parts[epsilon_idx + 1])
        except (ValueError, IndexError):
            print(f"⚠ Could not parse gamma/epsilon from {case_dir.name}")
            continue

        plot_case_comparison(case_dir, output_dir, gamma_val, epsilon_val, u_max, delta_x)

    print("\n" + "=" * 80)
    print("Reference-style analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
