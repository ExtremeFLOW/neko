#!/usr/bin/env python3
"""
Bound Violation Analysis for Advecting Drop Parameter Sweep

This script sweeps through all gamma_*_epsilon_* directories and computes
bound violations for the scalar field, which should be bounded in [0, 1].
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
except ImportError:
    print("Error: pysemtools not found. Please install pysemtools.")
    sys.exit(1)

comm = MPI.COMM_WORLD


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


def compute_bound_violation(field_data, scalar_min=0.0, scalar_max=1.0):
    """
    Compute bound violation statistics for a scalar field.

    Args:
        field_data: Field data array
        scalar_min: Minimum allowed value
        scalar_max: Maximum allowed value

    Returns:
        dict: Dictionary with violation statistics
    """
    all_scalar = field_data.flatten()

    # Calculate violation magnitudes
    upper_violation_mag = np.maximum(0, all_scalar - scalar_max)
    lower_violation_mag = np.maximum(0, scalar_min - all_scalar)
    total_violation_mag = upper_violation_mag + lower_violation_mag

    # Statistics
    max_violation = np.max(total_violation_mag)
    mean_violation = np.mean(total_violation_mag)

    # Counts
    total_points = len(all_scalar)
    upper_count = np.sum(upper_violation_mag > 0)
    lower_count = np.sum(lower_violation_mag > 0)
    violation_count = upper_count + lower_count
    violation_fraction = violation_count / total_points

    # Min/max values
    field_min = np.min(all_scalar)
    field_max = np.max(all_scalar)

    return {
        "max_violation": max_violation,
        "mean_violation": mean_violation,
        "violation_count": violation_count,
        "violation_fraction": violation_fraction,
        "upper_count": upper_count,
        "lower_count": lower_count,
        "field_min": field_min,
        "field_max": field_max,
        "total_points": total_points,
    }


def main():
    # Configuration
    base_dir = Path(
        "../data/advecting_drop_param_analysis_data"
    )  # Data directory containing gamma_* folders
    output_dir = Path(".")  # Current directory for outputs
    field_file_name = "field0.f00001"  # Final time step (t=1.0, with output_value=1.0)

    # Nondimensionalization constants
    u_max = 1.0  # Maximum velocity
    delta_x = 0.00357  # Average GLL spacing for 80x80 mesh on 2x2 domain

    print("=" * 80)
    print("Bound Violation Analysis for Advecting Drop Parameter Sweep")
    print("=" * 80)

    # Find all case directories
    case_dirs = sorted(base_dir.glob("gamma_*_epsilon_*"))

    if len(case_dirs) == 0:
        print("Error: No case directories found matching 'gamma_*_epsilon_*'")
        print(f"Searched in: {base_dir.resolve()}")
        sys.exit(1)

    print(f"\nFound {len(case_dirs)} case directories")

    # Storage for results
    results = []

    # Process each case
    for case_dir in case_dirs:
        dir_name = case_dir.name
        parsed = parse_directory_name(dir_name)

        if parsed is None:
            continue

        gamma_str, epsilon_str, gamma_float, epsilon_float = parsed

        field_file = case_dir / field_file_name

        if not field_file.exists():
            print(f"⚠ Warning: {field_file} not found, skipping...")
            continue

        try:
            # Read field data
            field_data = preadnek(str(field_file), comm)
            fld = field_c(comm, data=field_data)

            # Compute bound violations
            stats = compute_bound_violation(fld.fields["scal"][0])

            # Store results
            result = {
                "gamma": gamma_float,
                "epsilon": epsilon_float,
                "gamma_str": gamma_str,
                "epsilon_str": epsilon_str,
                "dir_name": dir_name,
                **stats,
            }
            results.append(result)

            # Print progress
            print(
                f"✓ {dir_name:40s} | γ={gamma_str:8s} ε={epsilon_str:8s} | "
                f"Max violation: {stats['max_violation']:.6f} | "
                f"Violations: {stats['violation_count']}/{stats['total_points']}"
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
    max_violations = np.array([r["max_violation"] for r in results])
    mean_violations = np.array([r["mean_violation"] for r in results])
    violation_fractions = np.array([r["violation_fraction"] for r in results])

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
        f"\nMax violation range: [{np.min(max_violations):.6f}, {np.max(max_violations):.6f}]"
    )

    # ========================================================================
    # PLOTTING
    # ========================================================================

    print("\nGenerating plots...")

    # Plot 1: Max violation vs gamma/u_max (for each epsilon/delta_x)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=150)

    for eps_nd in unique_epsilons_nd:
        mask = epsilons_nd == eps_nd
        gamma_subset = gammas_nd[mask]
        viol_subset = max_violations[mask]

        # Sort by gamma for proper line plotting
        sort_idx = np.argsort(gamma_subset)
        ax.loglog(
            gamma_subset[sort_idx],
            viol_subset[sort_idx],
            "o-",
            label=r"$\varepsilon/\Delta x$ = " + f"{eps_nd:.1f}",
            markersize=5,
            linewidth=1.5,
        )

    ax.set_xlabel(r"$\gamma / u_{\mathrm{max}}$", fontsize=12)
    ax.set_ylabel("Maximum Bound Violation", fontsize=12)
    ax.set_title(
        r"Maximum Scalar Field Bound Violation vs $\gamma/u_{\mathrm{max}}$ (for each $\varepsilon/\Delta x$)",
        fontsize=14,
    )
    ax.legend(loc="best", fontsize=9, ncol=2)
    ax.grid(True, which="both", ls="--", alpha=0.5)

    output_file = output_dir / "bound_violation_vs_gamma.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✓ Saved: {output_file}")
    plt.close()

    # Plot 2: Max violation vs epsilon/delta_x (for each gamma/u_max)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=150)

    for gam_nd in unique_gammas_nd:
        mask = gammas_nd == gam_nd
        eps_subset = epsilons_nd[mask]
        viol_subset = max_violations[mask]

        # Sort by epsilon
        sort_idx = np.argsort(eps_subset)
        ax.semilogy(
            eps_subset[sort_idx],
            viol_subset[sort_idx],
            "o-",
            label=r"$\gamma/u_{\mathrm{max}}$ = " + f"{gam_nd}",
            markersize=5,
            linewidth=1.5,
        )

    ax.set_xlabel(r"$\varepsilon / \Delta x$", fontsize=12)
    ax.set_ylabel("Maximum Bound Violation", fontsize=12)
    ax.set_title(
        r"Maximum Scalar Field Bound Violation vs $\varepsilon/\Delta x$ (for each $\gamma/u_{\mathrm{max}}$)",
        fontsize=14,
    )
    ax.legend(loc="best", fontsize=8, ncol=3)
    ax.grid(True, which="both", ls="--", alpha=0.5)

    output_file = output_dir / "bound_violation_vs_epsilon.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✓ Saved: {output_file}")
    plt.close()

    # Plot 3: 2D categorical heatmap of bound violations (nondimensional axes)
    # Filter out gamma=0 for log-scale plotting (log(0) is undefined)
    unique_gammas_nd_nonzero = unique_gammas_nd[unique_gammas_nd > 0]

    if len(unique_gammas_nd_nonzero) > 1 and len(unique_epsilons_nd) > 1:
        # Create 2D grid for violations
        violation_grid = np.zeros((len(unique_epsilons_nd), len(unique_gammas_nd_nonzero)))
        violation_grid[:] = np.nan

        # Create grid for stability status: 0=bounded, 1=unbounded, 2=unstable
        stability_grid = np.zeros((len(unique_epsilons_nd), len(unique_gammas_nd_nonzero)))
        stability_grid[:] = np.nan

        # Tolerance thresholds for classification
        tol_1e5 = 1e-5
        tol_1e3 = 1e-3
        unstable_threshold = 1.0

        for i, eps_nd in enumerate(unique_epsilons_nd):
            for j, gam_nd in enumerate(unique_gammas_nd_nonzero):
                mask = (epsilons_nd == eps_nd) & (gammas_nd == gam_nd)
                if np.any(mask):
                    idx = np.where(mask)[0][0]
                    violation = max_violations[idx]
                    violation_grid[i, j] = violation

                    # Classify stability with granular tolerance levels
                    if np.isnan(violation) or violation > unstable_threshold:
                        stability_grid[i, j] = 4  # Unstable
                    elif violation > tol_1e3:
                        stability_grid[i, j] = 3  # Unbounded (> 10^-3)
                    elif violation > tol_1e5:
                        stability_grid[i, j] = 2  # Within 10^-3
                    elif violation > 0:
                        stability_grid[i, j] = 1  # Within 10^-5
                    else:
                        stability_grid[i, j] = 0  # Bounded (no violation)

        # Count cases in each category
        n_bounded = np.sum(stability_grid == 0)
        n_within_1e5 = np.sum(stability_grid == 1)
        n_within_1e3 = np.sum(stability_grid == 2)
        n_unbounded = np.sum(stability_grid == 3)
        n_unstable = np.sum(stability_grid == 4)

        # Diagnostic output
        print(f"\n  Stability statistics for heatmap:")
        print(f"  Bounded (no violation): {n_bounded} cases")
        print(f"  Within 10^-5: {n_within_1e5} cases")
        print(f"  Within 10^-3: {n_within_1e3} cases")
        print(f"  Unbounded (> 10^-3): {n_unbounded} cases")
        print(f"  Unstable (blow-up): {n_unstable} cases")

        fig, ax = plt.subplots(1, 1, figsize=(12, 8), dpi=150)

        # Create discrete colormap with granular categories
        from matplotlib.colors import ListedColormap, BoundaryNorm

        # Black, Forest Green, Gold, Light Gray, Light Gray
        colors_list = ["#000000", "#228B22", "#FFD700", "#C0C0C0", "#C0C0C0"]
        cmap = ListedColormap(colors_list)
        boundaries = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
        norm = BoundaryNorm(boundaries, cmap.N)

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
                [
                    unique_epsilons_nd[0]
                    - (unique_epsilons_nd[1] - unique_epsilons_nd[0]) / 2
                ],
                (unique_epsilons_nd[:-1] + unique_epsilons_nd[1:]) / 2,
                [
                    unique_epsilons_nd[-1]
                    + (unique_epsilons_nd[-1] - unique_epsilons_nd[-2]) / 2
                ],
            ]
        )

        # Plot the base colors
        im = ax.pcolormesh(
            gamma_edges_nd,
            epsilon_edges_nd,
            stability_grid,
            cmap=cmap,
            norm=norm,
            shading="flat",
        )

        # Add hatching for unstable cases (denser pattern, more visible)
        for i, eps_nd in enumerate(unique_epsilons_nd):
            for j, gam_nd in enumerate(unique_gammas_nd_nonzero):
                if stability_grid[i, j] == 4:  # Unstable
                    # Draw hatching pattern
                    rect = plt.Rectangle(
                        (gamma_edges_nd[j], epsilon_edges_nd[i]),
                        gamma_edges_nd[j + 1] - gamma_edges_nd[j],
                        epsilon_edges_nd[i + 1] - epsilon_edges_nd[i],
                        fill=False,
                        hatch="////",
                        edgecolor="black",
                        linewidth=0,
                    )
                    ax.add_patch(rect)

        ax.set_xscale("log")
        ax.set_xlabel(r"$\gamma / u_{\mathrm{max}}$", fontsize=14, fontweight="bold")
        ax.set_ylabel(r"$\varepsilon / \Delta x$", fontsize=14, fontweight="bold")
        ax.set_title("Maximum Absolute Bound Violation", fontsize=14, fontweight="bold")

        # Add case count annotation
        total_cases = n_bounded + n_within_1e5 + n_within_1e3 + n_unbounded + n_unstable
        count_text = f"Total: {total_cases} cases"
        ax.text(
            0.02,
            0.02,
            count_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="gray"),
        )

        # Add legend instead of colorbar
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor="black", label="Bounded (no violation)"),
            Patch(facecolor="#228B22", label=r"Within $10^{-5}$"),
            Patch(facecolor="#FFD700", label=r"Within $10^{-3}$"),
            Patch(facecolor="#C0C0C0", label=r"Unbounded ($> 10^{-3}$)"),
            Patch(
                facecolor="#C0C0C0",
                hatch="////",
                edgecolor="black",
                label=f"Unstable (|viol| > {unstable_threshold})",
            ),
        ]
        ax.legend(handles=legend_elements, loc="upper right", fontsize=9)

        output_file = output_dir / "bound_violation_heatmap.png"
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"✓ Saved: {output_file}")
        plt.close()

    # Plot 4: Violation fraction vs gamma/u_max (for each epsilon/delta_x)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=150)

    for eps_nd in unique_epsilons_nd:
        mask = epsilons_nd == eps_nd
        gamma_subset = gammas_nd[mask]
        frac_subset = violation_fractions[mask] * 100  # Convert to percentage

        sort_idx = np.argsort(gamma_subset)
        ax.semilogx(
            gamma_subset[sort_idx],
            frac_subset[sort_idx],
            "o-",
            label=r"$\varepsilon/\Delta x$ = " + f"{eps_nd:.1f}",
            markersize=5,
            linewidth=1.5,
        )

    ax.set_xlabel(r"$\gamma / u_{\mathrm{max}}$", fontsize=12)
    ax.set_ylabel("Violation Fraction (%)", fontsize=12)
    ax.set_title(
        r"Percentage of Points Violating Bounds vs $\gamma/u_{\mathrm{max}}$",
        fontsize=14,
    )
    ax.legend(loc="best", fontsize=9, ncol=2)
    ax.grid(True, which="both", ls="--", alpha=0.5)

    output_file = output_dir / "violation_fraction_vs_gamma.png"
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

    # Find worst cases
    worst_idx = np.argmax(max_violations)
    best_idx = np.argmin(max_violations)

    print(f"\nWorst case (maximum violation):")
    print(f"  Directory: {results[worst_idx]['dir_name']}")
    print(f"  γ/u_max = {gammas_nd[worst_idx]}, ε/Δx = {epsilons_nd[worst_idx]:.2f}")
    print(f"  Max violation: {results[worst_idx]['max_violation']:.6f}")
    print(
        f"  Field range: [{results[worst_idx]['field_min']:.6f}, {results[worst_idx]['field_max']:.6f}]"
    )

    print(f"\nBest case (minimum violation):")
    print(f"  Directory: {results[best_idx]['dir_name']}")
    print(f"  γ/u_max = {gammas_nd[best_idx]}, ε/Δx = {epsilons_nd[best_idx]:.2f}")
    print(f"  Max violation: {results[best_idx]['max_violation']:.6f}")
    print(
        f"  Field range: [{results[best_idx]['field_min']:.6f}, {results[best_idx]['field_max']:.6f}]"
    )

    # Cases with no violations
    no_violation_cases = [
        (i, r) for i, r in enumerate(results) if r["max_violation"] < 1e-10
    ]
    if no_violation_cases:
        print(f"\nCases with no violations ({len(no_violation_cases)} total):")
        for i, r in no_violation_cases[:5]:  # Show first 5
            print(f"  γ/u_max = {gammas_nd[i]}, ε/Δx = {epsilons_nd[i]:.2f}")
        if len(no_violation_cases) > 5:
            print(f"  ... and {len(no_violation_cases) - 5} more")
    else:
        print("\nNo cases without violations found.")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
