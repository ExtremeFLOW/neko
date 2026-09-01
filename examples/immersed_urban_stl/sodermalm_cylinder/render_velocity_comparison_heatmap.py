#!/usr/bin/env python3
"""Render no-building/building terrain-relative velocity comparison panels."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle

from render_single_velocity_heatmap import fill_nearest, load_geojson, metadata_path, polygon_rings, smooth


def local_ring(ring: list[tuple[float, float]], center: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    pts = np.array([(x - center[0], y - center[1]) for x, y in ring], dtype=np.float64)
    return pts[:, 0], pts[:, 1]


def draw_context(ax, generated: Path, center: tuple[float, float], radius: float, *, buildings_alpha: float = 0.18) -> None:
    ax.add_patch(Circle((0.0, 0.0), radius, fill=False, edgecolor="#41545d", lw=0.65, alpha=0.45, zorder=3))
    land = generated / "sodermalm_osm_island_epsg3006.geojson"
    if land.exists():
        for ring in polygon_rings(load_geojson(land)["features"][0]["geometry"]):
            xs, ys = local_ring(ring, center)
            ax.plot(xs, ys, color="#112126", lw=0.85, alpha=0.9, zorder=4)
    buildings = generated / "sodermalm_buildings.geojson"
    if buildings.exists():
        for feature in load_geojson(buildings)["features"]:
            for ring in polygon_rings(feature["geometry"]):
                xs, ys = local_ring(ring, center)
                if len(xs) >= 3:
                    ax.fill(xs, ys, facecolor="#353535", edgecolor="none", alpha=buildings_alpha, zorder=5)


def load_speed(path: Path) -> tuple[np.lib.npyio.NpzFile, np.ndarray, np.ndarray]:
    data = np.load(path)
    speed_raw = data["speed"].astype(np.float32)
    disk = data["disk"].astype(bool)
    valid = data["valid"].astype(bool) & np.isfinite(speed_raw) & disk
    speed = smooth(fill_nearest(speed_raw, valid), disk, passes=2)
    speed = np.where(disk, speed, np.nan)
    return data, speed, valid


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit("usage: render_velocity_comparison_heatmap.py GENERATED NOBUILD_NPZ BUILD_NPZ OUT_PNG")
    generated = Path(sys.argv[1])
    no_data, no_speed, no_valid = load_speed(Path(sys.argv[2]))
    bld_data, bld_speed, bld_valid = load_speed(Path(sys.argv[3]))
    out = Path(sys.argv[4])

    meta = json.loads(metadata_path(generated).read_text(encoding="utf-8"))
    center = tuple(float(v) for v in meta["center_epsg3006"])
    radius = float(meta["radius_m"])
    disk = no_data["disk"].astype(bool) & bld_data["disk"].astype(bool)
    valid = no_valid & bld_valid & disk
    delta = bld_speed - no_speed
    extent = [
        float(np.nanmin(no_data["x"])),
        float(np.nanmax(no_data["x"])),
        float(np.nanmin(no_data["y"])),
        float(np.nanmax(no_data["y"])),
    ]

    speed_vals = np.concatenate([no_speed[valid], bld_speed[valid]])
    vmax = max(float(np.nanpercentile(speed_vals, 98.5)), 0.5)
    dmax = max(float(np.nanpercentile(np.abs(delta[valid]), 98.0)), 0.05)

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.9), dpi=180, constrained_layout=True)
    fig.patch.set_facecolor("white")
    speed_cmap = plt.get_cmap("turbo").copy()
    speed_cmap.set_bad("#d8edf6")
    diff_cmap = plt.get_cmap("RdBu_r").copy()
    diff_cmap.set_bad("#d8edf6")

    panels = [
        ("without buildings", no_speed, speed_cmap, dict(vmin=0.0, vmax=vmax)),
        ("with buildings", bld_speed, speed_cmap, dict(vmin=0.0, vmax=vmax)),
        ("with - without", delta, diff_cmap, dict(norm=TwoSlopeNorm(vmin=-dmax, vcenter=0.0, vmax=dmax))),
    ]
    images = []
    for ax, (title, values, cmap, kwargs) in zip(axes, panels):
        ax.set_facecolor("#d8edf6")
        img = ax.imshow(values, origin="upper", extent=extent, cmap=cmap, interpolation="bicubic", zorder=1, **kwargs)
        draw_context(ax, generated, center, radius, buildings_alpha=0.16 if "with" in title else 0.08)
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlim(-radius, radius)
        ax.set_ylim(-radius, radius)
        ax.set_xticks([])
        ax.set_yticks([])
        images.append(img)

    fig.colorbar(images[1], ax=axes[:2], shrink=0.82, pad=0.012, label="|u| (m/s)")
    fig.colorbar(images[2], ax=axes[2], shrink=0.82, pad=0.012, label="delta |u| (m/s)")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
