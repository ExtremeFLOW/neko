#!/usr/bin/env python3
"""Render terrain-relative velocity samples as a continuous top-down field."""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from sample_terrain_relative_velocity_npz import (
    color_velocity,
    load_geojson,
    metadata_path,
    polygon_boundary_rings,
    polygon_rings,
)


CANVAS_N = int(os.environ.get("CANVAS_N", "1600"))
BAR_W = 34
PAD = 34


def usage() -> None:
    raise SystemExit(
        "usage: render_velocity_field_heatmap.py GENERATED SAMPLE_NPZ OUT_PNG"
    )


def local_to_px(
    x: float,
    y: float,
    radius: float,
    offset_x: int = PAD,
    offset_y: int = PAD,
) -> tuple[int, int]:
    px = offset_x + int(round((x + radius) / (2.0 * radius) * (CANVAS_N - 1)))
    py = offset_y + int(round((radius - y) / (2.0 * radius) * (CANVAS_N - 1)))
    return px, py


def global_ring_to_px(
    ring: list[tuple[float, float]],
    center: tuple[float, float],
    radius: float,
) -> list[tuple[int, int]]:
    return [local_to_px(x - center[0], y - center[1], radius) for x, y in ring]


def draw_inlet_arc(
    draw: ImageDraw.ImageDraw,
    radius: float,
    wind_from_deg: float,
    arc_width_deg: float,
) -> None:
    pts = []
    for bearing in np.linspace(
        wind_from_deg - 0.5 * arc_width_deg,
        wind_from_deg + 0.5 * arc_width_deg,
        320,
    ):
        theta = math.radians(bearing)
        pts.append(local_to_px(radius * math.sin(theta), radius * math.cos(theta), radius))
    draw.line(pts, fill=(16, 94, 205, 255), width=max(6, CANVAS_N // 180), joint="curve")


def upscale_field(values: np.ndarray, disk: np.ndarray, valid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    finite = np.isfinite(values) & disk & valid
    filled = np.nan_to_num(values, nan=0.0).astype(np.float32)
    if finite.any():
        # Use nearest-neighbor growth only for gaps, then bicubic image resizing for display.
        known = finite.copy()
        for _ in range(max(values.shape)):
            if known.all():
                break
            acc = np.zeros_like(filled)
            count = np.zeros_like(filled)
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dx == 0 and dy == 0:
                        continue
                    k = np.roll(np.roll(known, dy, axis=0), dx, axis=1)
                    v = np.roll(np.roll(filled, dy, axis=0), dx, axis=1)
                    if dy < 0:
                        k[dy:, :] = False
                    elif dy > 0:
                        k[:dy, :] = False
                    if dx < 0:
                        k[:, dx:] = False
                    elif dx > 0:
                        k[:, :dx] = False
                    take = (~known) & k
                    acc[take] += v[take]
                    count[take] += 1.0
            grow = (~known) & (count > 0.0)
            filled[grow] = acc[grow] / count[grow]
            known[grow] = True

    src = Image.fromarray(filled.astype(np.float32), mode="F")
    smooth = np.asarray(src.resize((CANVAS_N, CANVAS_N), Image.Resampling.BICUBIC), dtype=np.float32)

    mask = Image.fromarray(disk.astype(np.uint8) * 255, mode="L")
    mask_big = np.asarray(mask.resize((CANVAS_N, CANVAS_N), Image.Resampling.BILINEAR), dtype=np.float32) / 255.0
    return smooth, mask_big


def add_colorbar(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    x0: int,
    y0: int,
    h: int,
    vmin: float,
    vmax: float,
) -> None:
    vals = np.linspace(vmax, vmin, h, dtype=np.float32)[:, None]
    bar = color_velocity(vals, vmin, vmax).reshape(h, 1, 3)
    bar = np.repeat(bar, BAR_W, axis=1)
    canvas.paste(Image.fromarray(bar, "RGB"), (x0, y0))
    draw.rectangle((x0, y0, x0 + BAR_W, y0 + h), outline=(20, 20, 20), width=1)
    for frac, label in ((0.0, f"{vmax:g}"), (0.5, f"{0.5 * (vmin + vmax):g}"), (1.0, f"{vmin:g}")):
        y = y0 + int(round(frac * h))
        draw.line((x0 + BAR_W, y, x0 + BAR_W + 8, y), fill=(20, 20, 20), width=1)
        draw.text((x0 + BAR_W + 12, y - 7), label, fill=(20, 20, 20))
    draw.text((x0, y0 - 22), "|u|", fill=(20, 20, 20))


def main() -> None:
    if len(sys.argv) != 4:
        usage()

    generated = Path(sys.argv[1])
    sample_path = Path(sys.argv[2])
    out = Path(sys.argv[3])
    meta = json.loads(metadata_path(generated).read_text())
    center = tuple(float(v) for v in meta["center_epsg3006"])
    radius = float(meta["radius_m"])
    wind_from = float(meta.get("inflow_from_degrees", meta.get("inflow_from_deg", 225.0)))
    arc_width = float(meta.get("inflow_arc_width_degrees", meta.get("inflow_arc_width_deg", 220.0)))

    sample = np.load(sample_path)
    speed = sample["speed"].astype(np.float32)
    disk = sample["disk"].astype(bool)
    valid = sample["valid"].astype(bool)
    field, mask = upscale_field(speed, disk, valid)

    vmin = float(os.environ.get("VMIN", "0.0"))
    vmax_override = os.environ.get("VMAX")
    visible_values = field[mask > 0.5]
    if visible_values.size == 0:
        raise RuntimeError("No visible velocity samples to render")
    vmax = float(vmax_override) if vmax_override else float(np.nanmax(visible_values))
    vmax = max(vmax, vmin + 1.0e-6)
    rgb = color_velocity(field, vmin, vmax)
    alpha = np.clip(mask * 255.0, 0, 255).astype(np.uint8)

    water = np.array([205, 226, 236], dtype=np.uint8)
    blended = (rgb.astype(np.float32) * (alpha[..., None] / 255.0) + water * (1.0 - alpha[..., None] / 255.0)).astype(np.uint8)

    canvas_w = CANVAS_N + 3 * PAD + BAR_W + 64
    canvas_h = CANVAS_N + 2 * PAD
    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
    canvas.paste(Image.fromarray(blended, "RGB"), (PAD, PAD))
    draw = ImageDraw.Draw(canvas, "RGBA")

    land_path = generated / "sodermalm_osm_island_epsg3006.geojson"
    if land_path.exists():
        land = load_geojson(land_path)["features"][0]["geometry"]
        for ring in polygon_boundary_rings(land):
            draw.line(global_ring_to_px(ring, center, radius), fill=(15, 25, 28, 215), width=max(2, CANVAS_N // 520), joint="curve")

    buildings_path = generated / "sodermalm_buildings.geojson"
    if buildings_path.exists() and os.environ.get("SHOW_BUILDINGS", "1") == "1":
        for feature in load_geojson(buildings_path)["features"]:
            for ring in polygon_rings(feature["geometry"]):
                pts = ring[:-1] if ring and ring[0] == ring[-1] else ring
                if len(pts) >= 3:
                    draw.polygon(
                        global_ring_to_px(pts, center, radius),
                        fill=(96, 96, 92, int(os.environ.get("BUILDING_ALPHA", "105"))),
                        outline=(45, 45, 42, 130),
                    )

    if os.environ.get("SHOW_INLET", "1") == "1":
        draw_inlet_arc(draw, radius, wind_from, arc_width)

    if os.environ.get("SHOW_COLORBAR", "1") == "1":
        add_colorbar(canvas, draw, CANVAS_N + 2 * PAD, PAD + CANVAS_N // 8, 3 * CANVAS_N // 4, vmin, vmax)

    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out)
    print(out)


if __name__ == "__main__":
    main()
