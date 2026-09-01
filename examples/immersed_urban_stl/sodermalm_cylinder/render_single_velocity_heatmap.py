#!/usr/bin/env python3
"""Render a single terrain-relative velocity heatmap for the urban gate cases."""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import PowerNorm
from matplotlib.patches import Circle


def load_geojson(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def metadata_path(generated: Path) -> Path:
    matches = sorted(p for p in generated.glob("*metadata.json") if not p.name.startswith("._"))
    if not matches:
        raise FileNotFoundError(f"No metadata JSON found in {generated}")
    return matches[0]


def polygon_rings(geometry: dict) -> list[list[tuple[float, float]]]:
    if geometry["type"] == "Polygon":
        polygons = [geometry["coordinates"]]
    elif geometry["type"] == "MultiPolygon":
        polygons = geometry["coordinates"]
    else:
        return []
    return [[(float(pt[0]), float(pt[1])) for pt in ring] for poly in polygons for ring in poly]


def fill_nearest(values: np.ndarray, valid: np.ndarray) -> np.ndarray:
    filled = np.nan_to_num(values, nan=0.0).astype(np.float32)
    known = valid.copy()
    for _ in range(max(values.shape)):
        if known.all():
            break
        acc = np.zeros_like(filled)
        count = np.zeros_like(filled)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dx == 0 and dy == 0:
                    continue
                shifted_known = np.roll(np.roll(known, dy, axis=0), dx, axis=1)
                shifted_vals = np.roll(np.roll(filled, dy, axis=0), dx, axis=1)
                if dy < 0:
                    shifted_known[dy:, :] = False
                elif dy > 0:
                    shifted_known[:dy, :] = False
                if dx < 0:
                    shifted_known[:, dx:] = False
                elif dx > 0:
                    shifted_known[:, :dx] = False
                take = (~known) & shifted_known
                acc[take] += shifted_vals[take]
                count[take] += 1.0
        grow = (~known) & (count > 0.0)
        filled[grow] = acc[grow] / count[grow]
        known[grow] = True
    return filled


def smooth(values: np.ndarray, valid: np.ndarray, passes: int = 4) -> np.ndarray:
    out = values.astype(np.float32, copy=True)
    for _ in range(passes):
        pad = np.pad(out, 1, mode="edge")
        candidate = (
            4.0 * pad[1:-1, 1:-1]
            + 2.0 * (pad[:-2, 1:-1] + pad[2:, 1:-1] + pad[1:-1, :-2] + pad[1:-1, 2:])
            + pad[:-2, :-2]
            + pad[:-2, 2:]
            + pad[2:, :-2]
            + pad[2:, 2:]
        ) / 16.0
        out[valid] = candidate[valid]
    return out


def local_ring(ring: list[tuple[float, float]], center: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    pts = np.array([(x - center[0], y - center[1]) for x, y in ring], dtype=np.float64)
    return pts[:, 0], pts[:, 1]


def draw_inlet_arc(ax, radius: float, wind_from_deg: float, arc_width_deg: float) -> None:
    angles = np.linspace(wind_from_deg - 0.5 * arc_width_deg, wind_from_deg + 0.5 * arc_width_deg, 240)
    theta = np.deg2rad(90.0 - angles)
    ax.plot(radius * np.cos(theta), radius * np.sin(theta), color="#1267d8", lw=5.0, solid_capstyle="round", zorder=8)

    flow_to = math.radians(90.0 - ((wind_from_deg + 180.0) % 360.0))
    dx = math.cos(flow_to)
    dy = math.sin(flow_to)
    ax.arrow(
        -0.33 * radius * dx,
        -0.33 * radius * dy,
        0.26 * radius * dx,
        0.26 * radius * dy,
        width=0.012 * radius,
        head_width=0.055 * radius,
        head_length=0.075 * radius,
        color="#10283a",
        alpha=0.88,
        length_includes_head=True,
        zorder=9,
    )


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: render_single_velocity_heatmap.py GENERATED SAMPLE_NPZ OUT_PNG")

    generated = Path(sys.argv[1])
    sample = np.load(sys.argv[2])
    out = Path(sys.argv[3])
    meta = json.loads(metadata_path(generated).read_text())
    radius = float(meta["radius_m"])
    center = tuple(float(v) for v in meta["center_epsg3006"])

    speed_raw = sample["speed"].astype(np.float32)
    disk = sample["disk"].astype(bool)
    valid = sample["valid"].astype(bool) & np.isfinite(speed_raw) & disk
    speed = smooth(fill_nearest(speed_raw, valid), disk, passes=3)
    speed = np.where(disk, speed, np.nan)

    finite = np.isfinite(speed) & disk
    vmax_override = os.environ.get("VMAX")
    vmax = float(vmax_override) if vmax_override else float(np.percentile(speed[finite], 99.2))
    vmax = max(vmax, 1.0)
    x = sample["x"]
    y = sample["y"]
    extent = [float(np.nanmin(x)), float(np.nanmax(x)), float(np.nanmin(y)), float(np.nanmax(y))]

    cmap = plt.get_cmap("turbo").copy()
    cmap.set_bad("#d8edf6")
    crop_half_width = float(os.environ.get("CROP_HALF_WIDTH", "0.0"))
    crop_x = float(os.environ.get("CROP_X", "0.0"))
    crop_y = float(os.environ.get("CROP_Y", "0.0"))
    show_inlet = os.environ.get("SHOW_INLET", "1") == "1"
    show_land = os.environ.get("SHOW_LAND", "1") == "1"
    show_buildings = os.environ.get("SHOW_BUILDINGS", "1") == "1"
    time_offset = float(os.environ.get("TIME_OFFSET", "0.0"))

    fig, ax = plt.subplots(figsize=(11.2, 10.0), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#d8edf6")
    image = ax.imshow(
        speed,
        origin="upper",
        extent=extent,
        cmap=cmap,
        norm=PowerNorm(gamma=0.58, vmin=0.0, vmax=vmax),
        interpolation="bicubic",
        zorder=1,
    )

    ax.add_patch(Circle((0.0, 0.0), radius, fill=False, edgecolor="#56656b", lw=0.8, alpha=0.45, zorder=3))
    land = generated / "sodermalm_osm_island_epsg3006.geojson"
    if show_land and land.exists():
        for ring in polygon_rings(load_geojson(land)["features"][0]["geometry"]):
            xs, ys = local_ring(ring, center)
            ax.plot(xs, ys, color="#112126", lw=1.2, alpha=0.82, zorder=4)

    buildings = generated / "sodermalm_buildings.geojson"
    if show_buildings and buildings.exists():
        for feature in load_geojson(buildings)["features"]:
            for ring in polygon_rings(feature["geometry"]):
                xs, ys = local_ring(ring, center)
                if len(xs) >= 3:
                    ax.fill(xs, ys, facecolor="#626262", edgecolor="#1d1d1d", lw=0.42, alpha=0.26, zorder=5)

    if show_inlet:
        wind_from = float(meta.get("inflow_from_degrees", meta.get("inflow_from_deg", 225.0)))
        arc_width = float(meta.get("inflow_arc_width_degrees", meta.get("inflow_arc_width_deg", 180.0)))
        draw_inlet_arc(ax, radius, wind_from, arc_width)
    ax.set_aspect("equal")
    if crop_half_width > 0.0:
        ax.set_xlim(crop_x - crop_half_width, crop_x + crop_half_width)
        ax.set_ylim(crop_y - crop_half_width, crop_y + crop_half_width)
    else:
        ax.set_xlim(-radius, radius)
        ax.set_ylim(-radius, radius)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Velocity magnitude at terrain + 10 m, t = {float(sample['time']) + time_offset:.1f} s")
    cbar = fig.colorbar(image, ax=ax, shrink=0.82, pad=0.025)
    cbar.set_label("|u| (m/s)")
    cbar.ax.tick_params(length=0)
    if show_inlet:
        ax.text(
            0.02,
            0.02,
            "blue arc: inlet, arrow: wind direction",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            color="#10283a",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 4},
            zorder=10,
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
