#!/usr/bin/env python3
"""
Create Animation of Advecting Drop

This script creates an animation showing the droplet advecting through
the domain and returning to its original position.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
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


def load_field_file(filepath, msh=None):
    """
    Load a field file and extract 2D slice at z=0.

    Args:
        filepath: Path to field file
        msh: Optional pre-loaded mesh (for speed)

    Returns:
        tuple: (time, x_2d, y_2d, phi_2d, msh)
    """
    data = preadnek(str(filepath), comm)

    # Load mesh if not provided
    if msh is None:
        msh = msh_c(comm, data=data)

    # Load field
    fld = field_c(comm, data=data)
    phi = fld.fields["scal"][0]
    time = data.time

    # Extract 2D slice at z=0
    z_tol = 0.1
    z_mask = np.abs(msh.z) < z_tol

    x_2d = msh.x[z_mask]
    y_2d = msh.y[z_mask]
    phi_2d = phi[z_mask]

    return time, x_2d, y_2d, phi_2d, msh


def create_frame_images(output_dir, frames_dir, mesh_size=None, gamma=None, epsilon=None):
    """
    Create individual frame images from simulation output.

    Args:
        output_dir: Directory containing field files
        frames_dir: Directory to save frame images
        mesh_size: Mesh size parameter (optional)
        gamma: Gamma parameter (optional)
        epsilon: Epsilon parameter (optional)
    """
    output_dir = Path(output_dir)
    frames_dir = Path(frames_dir)
    frames_dir.mkdir(exist_ok=True)

    # Find all field files
    field_files = sorted(output_dir.glob("field0.f*"))

    if len(field_files) == 0:
        print(f"Error: No field files found in {output_dir}")
        return

    print(f"Found {len(field_files)} field files")

    # Load mesh once
    print("Loading mesh...")
    data = preadnek(str(field_files[0]), comm)
    msh = msh_c(comm, data=data)
    print("Mesh loaded")

    # Process each frame
    print("\nCreating frames...")
    for i, field_file in enumerate(field_files):
        try:
            time, x_2d, y_2d, phi_2d, _ = load_field_file(field_file, msh)

            # Create figure
            fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=150)

            # Plot filled contours
            levels = np.linspace(0, 1, 21)  # 21 levels from 0 to 1
            tcf = ax.tricontourf(
                x_2d,
                y_2d,
                phi_2d,
                levels=levels,
                cmap="RdBu_r",
                vmin=0,
                vmax=1,
                extend="both",
            )

            # Plot interface contour (φ=0.5)
            ax.tricontour(
                x_2d, y_2d, phi_2d, levels=[0.5], colors="black", linewidths=2
            )

            ax.set_xlabel("x", fontsize=14)
            ax.set_ylabel("y", fontsize=14)
            ax.set_title(
                f"Advecting Drop - Time: {time:.3f}", fontsize=16, fontweight="bold"
            )
            ax.set_aspect("equal")
            ax.set_xlim([0, 2])
            ax.set_ylim([0, 2])

            # Add colorbar with Greek phi symbol
            cbar = plt.colorbar(
                tcf, ax=ax, label="φ", ticks=np.linspace(0, 1, 11)
            )

            # Add grid
            
            plt.tight_layout()

            # Add parameter info below the plot (can be cropped out)
            param_text = []
            if mesh_size is not None:
                param_text.append(f"mesh: {mesh_size}")
            if gamma is not None:
                param_text.append(f"γ: {gamma}")
            if epsilon is not None:
                param_text.append(f"ε: {epsilon}")

            if param_text:
                fig.text(
                    0.5, 0.02, "  |  ".join(param_text),
                    ha="center", va="bottom", fontsize=10,
                    bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.8)
                )

            # Save frame
            frame_file = frames_dir / f"frame_{i:04d}.png"
            plt.savefig(frame_file, dpi=150, bbox_inches="tight")
            plt.close()

            print(f"  Frame {i+1}/{len(field_files)}: t={time:.3f} → {frame_file.name}")

        except Exception as e:
            print(f"  Error processing {field_file.name}: {e}")

    print(f"\n✓ All frames saved to: {frames_dir.resolve()}")


def create_animation(output_dir, animation_file, mesh_size=None, gamma=None, epsilon=None):
    """
    Create animated video from simulation output.

    Args:
        output_dir: Directory containing field files
        animation_file: Output animation filename (e.g., 'advecting_drop.mp4')
        mesh_size: Mesh size parameter (optional)
        gamma: Gamma parameter (optional)
        epsilon: Epsilon parameter (optional)
    """
    output_dir = Path(output_dir)

    # Find all field files
    field_files = sorted(output_dir.glob("field0.f*"))

    if len(field_files) == 0:
        print(f"Error: No field files found in {output_dir}")
        return

    print(f"Found {len(field_files)} field files")

    # Load mesh once
    print("Loading mesh...")
    data = preadnek(str(field_files[0]), comm)
    msh = msh_c(comm, data=data)
    print("Mesh loaded")

    # Load all frames
    print("\nLoading all frames...")
    frames_data = []
    for i, field_file in enumerate(field_files):
        try:
            time, x_2d, y_2d, phi_2d, _ = load_field_file(field_file, msh)
            frames_data.append((time, x_2d, y_2d, phi_2d))
            print(f"  Loaded {i+1}/{len(field_files)}: t={time:.3f}")
        except Exception as e:
            print(f"  Error loading {field_file.name}: {e}")

    if len(frames_data) == 0:
        print("Error: No frames loaded successfully")
        return

    print(f"\n✓ Loaded {len(frames_data)} frames")

    # Create animation
    print("\nCreating animation...")

    # Build parameter text
    param_text = []
    if mesh_size is not None:
        param_text.append(f"mesh: {mesh_size}")
    if gamma is not None:
        param_text.append(f"γ: {gamma}")
    if epsilon is not None:
        param_text.append(f"ε: {epsilon}")
    param_str = "  |  ".join(param_text) if param_text else None

    fig, ax = plt.subplots(1, 1, figsize=(10, 8), dpi=150)

    # Create a persistent colorbar using ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    norm = Normalize(vmin=0, vmax=1)
    sm = ScalarMappable(cmap="RdBu_r", norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label="φ", ticks=np.linspace(0, 1, 11))

    # Add parameter info below the plot (persistent, can be cropped out)
    if param_str:
        fig.text(
            0.5, 0.02, param_str,
            ha="center", va="bottom", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.8)
        )

    def animate(frame_idx):
        """Animation update function."""
        ax.clear()

        time, x_2d, y_2d, phi_2d = frames_data[frame_idx]

        # Plot filled contours
        levels = np.linspace(0, 1, 21)  # 21 levels from 0 to 1
        tcf = ax.tricontourf(
            x_2d,
            y_2d,
            phi_2d,
            levels=levels,
            cmap="RdBu_r",
            vmin=0,
            vmax=1,
            extend="both",
        )

        # Plot interface contour (φ=0.5)
        ax.tricontour(x_2d, y_2d, phi_2d, levels=[0.5], colors="black", linewidths=2)

        ax.set_xlabel("x", fontsize=14)
        ax.set_ylabel("y", fontsize=14)
        ax.set_title(
            f"Advecting Drop - Time: {time:.3f}", fontsize=16, fontweight="bold"
        )
        ax.set_aspect("equal")
        ax.set_xlim([0, 2])
        ax.set_ylim([0, 2])
        
        return (ax,)

    # Create animation
    anim = FuncAnimation(
        fig, animate, frames=len(frames_data), interval=200, blit=False, repeat=True
    )

    # Save animation
    try:
        writer = FFMpegWriter(fps=5, bitrate=2000)
        anim.save(animation_file, writer=writer)
        print(f"\n✓ Animation saved to: {Path(animation_file).resolve()}")
    except Exception as e:
        print(f"\nError saving animation: {e}")
        print("Note: ffmpeg must be installed to save animations")
        print("Try: brew install ffmpeg (macOS) or apt-get install ffmpeg (Linux)")

    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Create animation of advecting drop simulation"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="../visualization_output",
        help="Directory containing simulation output files",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["frames", "animation"],
        default="frames",
        help="Output mode: frames (individual images) or animation (video)",
    )
    parser.add_argument(
        "--frames-dir",
        type=str,
        default="animation_frames",
        help="Directory to save individual frames (for frames mode)",
    )
    parser.add_argument(
        "--animation-file",
        type=str,
        default="advecting_drop.mp4",
        help="Output animation filename (for animation mode)",
    )
    parser.add_argument(
        "--mesh-size",
        type=str,
        default=None,
        help="Mesh size parameter to display on frames",
    )
    parser.add_argument(
        "--gamma",
        type=str,
        default=None,
        help="Gamma parameter to display on frames",
    )
    parser.add_argument(
        "--epsilon",
        type=str,
        default=None,
        help="Epsilon parameter to display on frames",
    )

    args = parser.parse_args()

    print("=" * 80)
    print("Advecting Drop Animation Creator")
    print("=" * 80)
    print(f"Mode: {args.mode}")
    print(f"Output directory: {Path(args.output_dir).resolve()}")
    print("=" * 80)

    if args.mode == "frames":
        create_frame_images(
            args.output_dir, args.frames_dir,
            mesh_size=args.mesh_size, gamma=args.gamma, epsilon=args.epsilon
        )
    else:
        create_animation(
            args.output_dir, args.animation_file,
            mesh_size=args.mesh_size, gamma=args.gamma, epsilon=args.epsilon
        )

    print("\n" + "=" * 80)
    print("Done!")
    print("=" * 80)


if __name__ == "__main__":
    main()
