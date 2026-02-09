#!/usr/bin/env python3
"""
Shape Error Analysis for Zalesak Disk Parameter Sweep

This script sweeps through all gamma_*_epsilon_* directories and computes
shape error, which measures how much the solution differs from the initial condition.
The error is the L1 norm of the difference, integrated over the domain.
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
    from pysemtools.datatypes.coef import Coef as coef_c
except ImportError:
    print("Error: pysemtools not found. Please install pysemtools.")
    sys.exit(1)

comm = MPI.COMM_WORLD


def initial_condition(x, y, eps):
    """
    Generate the initial condition for Zalesak disk with notch.

    Args:
        x: x-coordinates
        y: y-coordinates
        eps: interface thickness parameter

    Returns:
        field: Initial condition field
    """
    # Disk parameters
    center_x = 0.5
    center_y = 0.75
    radius = 0.15
    notch_width = 0.05
    notch_height = 0.05

    # Create the disk with tanh profile
    rad = np.sqrt((x - center_x) ** 2 + (y - center_y) ** 2)
    field = 0.5 * (1 + np.tanh((radius - rad) / (2 * eps)))

    # Cut out the notch with smooth tanh transition
    notch_x = 0.5 * (1 + np.tanh((np.abs(x - center_x) - notch_width / 2) / (2 * eps)))
    notch_y = 0.5 * (1 + np.tanh((y - (center_y + notch_height)) / (2 * eps)))
    field = field * (notch_x + (1 - notch_x) * notch_y)

    return field


def parse_directory_name(dir_name):
    """
    Parse directory name like 'gamma_0.01_epsilon_0.02' to extract parameters.

    Returns:
        tuple: (gamma_str, epsilon_str, gamma_float, epsilon_float)
    """
    parts = dir_name.split("_")
    try:
        # Find gamma and epsilon values
        gamma_idx = parts.index("gamma")
        epsilon_idx = parts.index("epsilon")

        gamma_str = parts[gamma_idx + 1]
        epsilon_str = parts[epsilon_idx + 1]

        gamma_float = float(gamma_str)
        epsilon_float = float(epsilon_str)

        return gamma_str, epsilon_str, gamma_float, epsilon_float
    except (ValueError, IndexError) as e:
        print(f"Warning: Could not parse directory name '{dir_name}': {e}")
        return None


def compute_shape_error(initial_field, final_field, coef):
    """
    Compute shape error as L1 norm of difference between final and initial fields.

    Args:
        initial_field: Initial condition field
        final_field: Final field from simulation
        coef: Coefficient object containing mass matrix (B)

    Returns:
        float: Integrated shape error
    """
    # Compute absolute difference
    diff = np.abs(final_field - initial_field)

    # Integrate using mass matrix (L1 norm)
    error = np.sum(diff * coef.B)

    return error


def main():
    # Configuration
    base_dir = Path(
        "../data/zalesak_param_analysis_data"
    )  # Data directory containing gamma_* folders
    output_dir = Path(".")  # Current directory for outputs
    field_file_name = "field0.f00002"  # Final time step
    initial_file_name = "field0.f00000"  # Initial time step

    # Nondimensionalization constants
    u_max = 1.0  # Maximum velocity
    delta_x = 0.00357  # Average GLL spacing for 40x40 mesh

    print("=" * 80)
    print("Shape Error Analysis for Zalesak Disk Parameter Sweep")
    print("=" * 80)

    # Find all case directories
    case_dirs = sorted(base_dir.glob("gamma_*_epsilon_*"))

    if len(case_dirs) == 0:
        print("Error: No case directories found matching 'gamma_*_epsilon_*'")
        print(f"Searched in: {base_dir.resolve()}")
        sys.exit(1)

    print(f"\nFound {len(case_dirs)} case directories")

    # We need to load mesh and coefficients once (assuming same mesh for all cases)
    # Load from first available case
    first_case = case_dirs[0]
    first_initial_file = first_case / initial_file_name

    if not first_initial_file.exists():
        print(f"Error: Could not find initial file in {first_case}")
        sys.exit(1)

    print(f"\nLoading mesh and coefficients from {first_case.name}...")
    xyz_info = preadnek(str(first_initial_file), comm)
    msh = msh_c(comm, data=xyz_info)
    coef = coef_c(msh, comm)
    print("✓ Mesh loaded")

    # Storage for results
    results = []

    # Process each case
    print("\nProcessing cases...")
    for case_dir in case_dirs:
        dir_name = case_dir.name
        parsed = parse_directory_name(dir_name)

        if parsed is None:
            continue

        gamma_str, epsilon_str, gamma_float, epsilon_float = parsed

        initial_file = case_dir / initial_file_name
        final_file = case_dir / field_file_name

        if not initial_file.exists():
            print(f"⚠ Warning: {initial_file} not found, skipping...")
            continue

        if not final_file.exists():
            print(f"⚠ Warning: {final_file} not found, skipping...")
            continue

        try:
            # Load initial field to get mesh coordinates (in case they differ slightly)
            initial_data = preadnek(str(initial_file), comm)
            msh_case = msh_c(comm, data=initial_data)

            # Generate analytical initial condition using epsilon from this case
            ic = initial_condition(msh_case.x, msh_case.y, epsilon_float)

            # Load final field
            final_data = preadnek(str(final_file), comm)
            fld = field_c(comm, data=final_data)

            # Compute shape error
            shape_error = compute_shape_error(ic, fld.fields["scal"][0], coef)

            # Store results
            result = {
                "gamma": gamma_float,
                "epsilon": epsilon_float,
                "gamma_str": gamma_str,
                "epsilon_str": epsilon_str,
                "dir_name": dir_name,
                "shape_error": shape_error,
            }
            results.append(result)

            # Print progress
            print(
                f"✓ {dir_name:40s} | γ={gamma_str:8s} ε={epsilon_str:8s} | "
                f"Shape error: {shape_error:.6e}"
            )

        except Exception as e:
            print(f"✗ Error processing {dir_name}: {e}")
            continue

    if len(results) == 0:
        print("\nError: No results collected. Check that field files exist.")
        sys.exit(1)

    print(f"\n{'=' * 80}")
    print(f"Successfully processed {len(results)} cases")
    print(f"{'=' * 80}\n")

    # Convert to arrays for analysis
    gammas = np.array([r["gamma"] for r in results])
    epsilons = np.array([r["epsilon"] for r in results])
    shape_errors = np.array([r["shape_error"] for r in results])

    # Nondimensional arrays
    gammas_nd = gammas / u_max  # γ/u_max (same values since u_max=1)
    epsilons_nd = epsilons / delta_x  # ε/Δx

    # Get unique gamma and epsilon values (nondimensional)
    unique_gammas_nd = np.unique(gammas_nd)
    unique_epsilons_nd = np.unique(epsilons_nd)

    print(f"Unique γ/u_max values: {len(unique_gammas_nd)}")
    print(f"Unique ε/Δx values: {len(unique_epsilons_nd)}")
    print(f"\nγ/u_max range: [{np.min(gammas_nd)}, {np.max(gammas_nd)}]")
    print(f"ε/Δx range: [{np.min(epsilons_nd):.2f}, {np.max(epsilons_nd):.2f}]")
    print(
        f"\nShape error range: [{np.min(shape_errors):.6e}, {np.max(shape_errors):.6e}]"
    )

    # ========================================================================
    # PLOTTING
    # ========================================================================

    print("\nGenerating plots...")

    # Plot 1: Shape error vs gamma/u_max (for each epsilon/delta_x)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=150)

    for eps_nd in unique_epsilons_nd:
        mask = epsilons_nd == eps_nd
        gamma_subset = gammas_nd[mask]
        error_subset = shape_errors[mask]

        # Sort by gamma for proper line plotting
        sort_idx = np.argsort(gamma_subset)
        ax.loglog(
            gamma_subset[sort_idx],
            error_subset[sort_idx],
            "o-",
            label=r"$\varepsilon/\Delta x$ = " + f"{eps_nd:.1f}",
            markersize=5,
            linewidth=1.5,
        )

    ax.set_xlabel(r"$\gamma / u_{\mathrm{max}}$", fontsize=12)
    ax.set_ylabel("Shape Error (L1 norm)", fontsize=12)
    ax.set_title(r"Shape Error vs $\gamma/u_{\mathrm{max}}$ (for each $\varepsilon/\Delta x$)", fontsize=14)
    ax.legend(loc="best", fontsize=9, ncol=2)
    ax.grid(True, which="both", ls="--", alpha=0.5)

    output_file = output_dir / "shape_error_vs_gamma.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✓ Saved: {output_file}")
    plt.close()

    # Plot 2: Shape error vs epsilon/delta_x (for each gamma/u_max)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=150)

    for gam_nd in unique_gammas_nd:
        mask = gammas_nd == gam_nd
        eps_subset = epsilons_nd[mask]
        error_subset = shape_errors[mask]

        # Sort by epsilon
        sort_idx = np.argsort(eps_subset)
        ax.semilogy(
            eps_subset[sort_idx],
            error_subset[sort_idx],
            "o-",
            label=r"$\gamma/u_{\mathrm{max}}$ = " + f"{gam_nd}",
            markersize=5,
            linewidth=1.5,
        )

    ax.set_xlabel(r"$\varepsilon / \Delta x$", fontsize=12)
    ax.set_ylabel("Shape Error (L1 norm)", fontsize=12)
    ax.set_title(r"Shape Error vs $\varepsilon/\Delta x$ (for each $\gamma/u_{\mathrm{max}}$)", fontsize=14)
    ax.legend(loc="best", fontsize=8, ncol=3)
    ax.grid(True, which="both", ls="--", alpha=0.5)

    output_file = output_dir / "shape_error_vs_epsilon.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✓ Saved: {output_file}")
    plt.close()

    # Plot 3: 2D heatmap of shape errors (nondimensional axes)
    # Filter out gamma=0 for log-scale plotting (log(0) is undefined)
    unique_gammas_nd_nonzero = unique_gammas_nd[unique_gammas_nd > 0]

    if len(unique_gammas_nd_nonzero) > 1 and len(unique_epsilons_nd) > 1:
        # Create 2D grid
        error_grid = np.zeros((len(unique_epsilons_nd), len(unique_gammas_nd_nonzero)))
        error_grid[:] = np.nan

        for i, eps_nd in enumerate(unique_epsilons_nd):
            for j, gam_nd in enumerate(unique_gammas_nd_nonzero):
                mask = (epsilons_nd == eps_nd) & (gammas_nd == gam_nd)
                if np.any(mask):
                    error_grid[i, j] = shape_errors[mask][0]

        fig, ax = plt.subplots(1, 1, figsize=(12, 8), dpi=150)

        # Create custom colormap: grayscale from white (low error) to black (high error)
        from matplotlib.colors import ListedColormap, LogNorm

        # Find min and max errors for scaling
        min_error = np.nanmin(error_grid)
        max_error = np.nanmax(error_grid)

        # Create grayscale colormap
        # Black for low error, white for high error
        colors_list = []

        # Add gradient from black to white for increasing errors
        n_colors = 100
        for i in range(n_colors):
            intensity = i / (n_colors - 1)  # 0 to 1
            # Grayscale: from black (0, 0, 0) to white (255, 255, 255)
            gray_value = int(255 * intensity)
            colors_list.append(f"#{gray_value:02x}{gray_value:02x}{gray_value:02x}")

        cmap = ListedColormap(colors_list)

        # Use logarithmic normalization for colorbar
        # Set colorbar limits to span the full data range
        colorbar_min = 5e-5  # Lower limit (just below data min ~6e-5)
        colorbar_max = 1.5e-2  # Upper limit (just above data max ~1.1e-2)

        # Use the actual data range or the specified limits, whichever is more restrictive
        vmin_use = max(min_error, colorbar_min) if min_error > 0 else colorbar_min
        vmax_use = colorbar_max

        # Use logarithmic normalization
        if max_error > min_error and min_error > 0:
            norm = LogNorm(vmin=vmin_use, vmax=vmax_use)
        else:
            # Fallback for edge cases
            norm = None

        # Print diagnostic information
        n_below = np.sum(error_grid < vmin_use)
        n_above = np.sum(error_grid > vmax_use)
        n_within = np.sum((error_grid >= vmin_use) & (error_grid <= vmax_use))
        n_total = n_below + n_above + n_within

        print(f"  Data range: [{min_error:.3e}, {max_error:.3e}]")
        print(f"  Colorbar range: [{vmin_use:.3e}, {vmax_use:.3e}]")
        print(
            f"  Cases: {n_below} below, {n_within} within, {n_above} above colorbar range (total: {n_total})"
        )

        # Use pcolormesh for better control
        # Create mesh grid for gamma/u_max and epsilon/delta_x (nondimensional)
        gamma_edges_nd = np.concatenate(
            [
                [unique_gammas_nd_nonzero[0] * 0.9],  # Extend slightly beyond
                (unique_gammas_nd_nonzero[:-1] + unique_gammas_nd_nonzero[1:]) / 2,  # Midpoints
                [unique_gammas_nd_nonzero[-1] * 1.1],
            ]
        )
        epsilon_edges_nd = np.concatenate(
            [
                [unique_epsilons_nd[0] - (unique_epsilons_nd[1] - unique_epsilons_nd[0]) / 2],
                (unique_epsilons_nd[:-1] + unique_epsilons_nd[1:]) / 2,
                [unique_epsilons_nd[-1] + (unique_epsilons_nd[-1] - unique_epsilons_nd[-2]) / 2],
            ]
        )

        im = ax.pcolormesh(
            gamma_edges_nd,
            epsilon_edges_nd,
            error_grid,
            cmap=cmap,
            norm=norm,
            shading="flat",
        )

        ax.set_xscale("log")
        ax.set_xlabel(r"$\gamma / u_{\mathrm{max}}$", fontsize=14, fontweight="bold")
        ax.set_ylabel(r"$\varepsilon / \Delta x$", fontsize=14, fontweight="bold")
        ax.set_title(
            "Shape Error Map: Black=Low Error, White=High Error",
            fontsize=14,
            fontweight="bold",
        )

        # Add case count annotation
        count_text = f"Total: {n_total} cases\n{n_below} below | {n_within} within | {n_above} above colorbar range"
        ax.text(
            0.02,
            0.02,
            count_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="gray"),
        )

        cbar = plt.colorbar(im, ax=ax, extend="both")
        cbar.set_label("Shape Error (L1 norm) - Log Scale", fontsize=12)

        output_file = output_dir / "shape_error_heatmap.png"
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"✓ Saved: {output_file}")
        plt.close()

    # ========================================================================
    # SUMMARY STATISTICS
    # ========================================================================

    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    # Find worst and best cases
    worst_idx = np.argmax(shape_errors)
    best_idx = np.argmin(shape_errors)

    print(f"\nWorst case (maximum shape error):")
    print(f"  Directory: {results[worst_idx]['dir_name']}")
    print(f"  γ/u_max = {results[worst_idx]['gamma']/u_max}, ε/Δx = {results[worst_idx]['epsilon']/delta_x:.2f}")
    print(f"  Shape error: {results[worst_idx]['shape_error']:.6e}")

    print(f"\nBest case (minimum shape error):")
    print(f"  Directory: {results[best_idx]['dir_name']}")
    print(f"  γ/u_max = {results[best_idx]['gamma']/u_max}, ε/Δx = {results[best_idx]['epsilon']/delta_x:.2f}")
    print(f"  Shape error: {results[best_idx]['shape_error']:.6e}")

    # Statistical summary
    print(f"\nError statistics:")
    print(f"  Mean: {np.mean(shape_errors):.6e}")
    print(f"  Median: {np.median(shape_errors):.6e}")
    print(f"  Std dev: {np.std(shape_errors):.6e}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
