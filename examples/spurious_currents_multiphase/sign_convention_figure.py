"""
Generate a diagram explaining the curvature sign convention bug
in the CSF surface tension model.

The key issue:
  - Normal n_hat = grad(phi)/|grad(phi)| points INWARD (toward phi=1 region)
  - div(n_hat) for an inward-pointing normal on a circle is NEGATIVE: -1/R
  - The surface tension force F = sigma * kappa * grad(phi)
  - If kappa = div(n_hat) = -1/R < 0, then F points outward (WRONG for surface tension)
  - Fix: kappa = -div(n_hat) = +1/R > 0, so F points inward (CORRECT)

Usage:
    MPLBACKEND=Agg python sign_convention_figure.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

# ── Parameters ──────────────────────────────────────────────────────────────
cx, cy = 0.5, 0.5   # circle centre
R = 0.20             # circle radius
n_arrows = 8         # number of arrows around circle


def draw_panel(ax, title, force_color, force_label_suffix,
               kappa_text, force_text, force_direction):
    """
    Draw one panel of the sign-convention diagram.

    Parameters
    ----------
    force_direction : str
        'outward' or 'inward' -- direction of the surface tension force arrows.
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)

    # ── Background labels ───────────────────────────────────────────────
    ax.text(0.08, 0.92, r'$\varphi = 0$', fontsize=13, color='0.35',
            ha='left', va='top')

    # ── Drop (filled circle) ────────────────────────────────────────────
    circle_fill = plt.Circle((cx, cy), R, fc='#d0e8ff', ec='#2277bb',
                              lw=2.0, zorder=2)
    ax.add_patch(circle_fill)
    ax.text(cx, cy, r'$\varphi = 1$', fontsize=12, ha='center', va='center',
            color='#1155aa', zorder=3)

    # ── Arrows ──────────────────────────────────────────────────────────
    angles = np.linspace(0, 2 * np.pi, n_arrows, endpoint=False)
    arrow_len_normal = 0.09    # length of n-hat arrows
    arrow_len_force  = 0.10    # length of force arrows

    for i, theta in enumerate(angles):
        # Point on the circle
        px = cx + R * np.cos(theta)
        py = cy + R * np.sin(theta)

        # Unit outward radial direction (away from centre)
        ex = np.cos(theta)
        ey = np.sin(theta)

        # n-hat points INWARD (toward phi=1, i.e. toward centre)
        nx = -ex
        ny = -ey

        # Draw n-hat arrow (black, thin) -- starts at circle, points inward
        ax.annotate('',
                    xy=(px + arrow_len_normal * nx,
                        py + arrow_len_normal * ny),
                    xytext=(px, py),
                    arrowprops=dict(arrowstyle='->', color='black',
                                    lw=1.4, shrinkA=0, shrinkB=0),
                    zorder=4)

        # Force arrow direction
        if force_direction == 'outward':
            fx, fy = ex, ey          # outward (wrong)
        else:
            fx, fy = -ex, -ey        # inward  (correct)

        # Draw force arrow (coloured, thicker)
        # For outward: start just outside circle, point outward
        # For inward: start well outside circle, point inward toward circle
        #   -- so they don't overlap with the n-hat arrows
        if force_direction == 'outward':
            offset = 0.015
            start_x = px + offset * ex
            start_y = py + offset * ey
            end_x = start_x + arrow_len_force * fx
            end_y = start_y + arrow_len_force * fy
        else:
            # Start from outside, arrow tip near the circle surface
            gap = 0.02
            start_x = px + (arrow_len_force + gap) * ex
            start_y = py + (arrow_len_force + gap) * ey
            end_x = px + gap * ex
            end_y = py + gap * ey

        ax.annotate('',
                    xy=(end_x, end_y),
                    xytext=(start_x, start_y),
                    arrowprops=dict(arrowstyle='->', color=force_color,
                                    lw=2.2, shrinkA=0, shrinkB=0),
                    zorder=5)

    # ── Arrow legends (manual) ──────────────────────────────────────────
    # n-hat label -- attach to the arrow near theta=pi/4
    ref_theta = angles[1]  # second arrow
    ref_px = cx + R * np.cos(ref_theta)
    ref_py = cy + R * np.sin(ref_theta)
    label_x = ref_px - arrow_len_normal * np.cos(ref_theta) - 0.02
    label_y = ref_py - arrow_len_normal * np.sin(ref_theta) + 0.04
    ax.text(label_x, label_y,
            r'$\hat{n} = \frac{\nabla\varphi}{|\nabla\varphi|}$',
            fontsize=10, color='black', ha='center', va='bottom',
            zorder=6,
            bbox=dict(fc='white', ec='none', alpha=0.85, pad=1))

    # F_ST label -- place near the bottom of the circle
    # For outward: label goes outside, beyond the outward force arrow tip
    # For inward: label goes outside too (between circle and axis edge) for clarity
    ref_theta2 = angles[n_arrows // 2]  # bottom arrow (theta ~ pi, pointing down-ish)
    ref_px2 = cx + R * np.cos(ref_theta2)
    ref_py2 = cy + R * np.sin(ref_theta2)
    if force_direction == 'outward':
        fl_x = ref_px2 + (arrow_len_force + 0.05) * np.cos(ref_theta2)
        fl_y = ref_py2 + (arrow_len_force + 0.05) * np.sin(ref_theta2)
    else:
        # Place label outside the circle even for inward arrows
        fl_x = ref_px2 + 0.14 * np.cos(ref_theta2)
        fl_y = ref_py2 + 0.14 * np.sin(ref_theta2)
    ax.text(fl_x, fl_y,
            r'$\mathbf{F}_{\mathrm{ST}}$',
            fontsize=11, fontweight='bold', color=force_color,
            ha='center', va='center', zorder=6,
            bbox=dict(fc='white', ec='none', alpha=0.85, pad=1))

    # ── Equations at bottom ─────────────────────────────────────────────
    ax.text(0.5, 0.09, kappa_text,
            fontsize=11, ha='center', va='center',
            transform=ax.transAxes,
            bbox=dict(fc='#f5f5f5', ec='0.7', boxstyle='round,pad=0.4'))
    ax.text(0.5, 0.02, force_text,
            fontsize=11, ha='center', va='center',
            transform=ax.transAxes, color=force_color, fontweight='bold')


# ── Create figure ───────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
fig.patch.set_facecolor('white')

# Left panel -- the bug
draw_panel(
    axes[0],
    title=r'Bug:  $\kappa = \nabla\!\cdot\!\hat{n} = -5$',
    force_color='#cc2222',
    force_label_suffix='(wrong)',
    kappa_text=r'$\kappa = \nabla\!\cdot\!\hat{n} = -\frac{1}{R} = -5.0$',
    force_text='$\\mathbf{F} = \\sigma\\,\\kappa\\,\\nabla\\varphi \\;\\;\\rightarrow\\;\\;$ outward  \u2717',
    force_direction='outward',
)

# Right panel -- the fix
draw_panel(
    axes[1],
    title=r'Fix:  $\kappa = -\nabla\!\cdot\!\hat{n} = +5$',
    force_color='#22aa22',
    force_label_suffix='(correct)',
    kappa_text=r'$\kappa = -\nabla\!\cdot\!\hat{n} = +\frac{1}{R} = +5.0$',
    force_text='$\\mathbf{F} = \\sigma\\,\\kappa\\,\\nabla\\varphi \\;\\;\\rightarrow\\;\\;$ inward  \u2713',
    force_direction='inward',
)

fig.tight_layout(pad=2.0)

# ── Save ────────────────────────────────────────────────────────────────────
import os
out_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(out_dir, 'sign_convention_diagram.png')
fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved: {out_path}")
