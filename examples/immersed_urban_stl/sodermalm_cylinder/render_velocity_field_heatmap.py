#!/usr/bin/env python3
"""Render terrain-relative velocity samples as a continuous top-down field."""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from matplotlib import colormaps, font_manager
from PIL import Image, ImageChops, ImageDraw, ImageFont

from sample_terrain_relative_velocity_npz import (
    load_geojson,
    metadata_path,
    polygon_boundary_rings,
)


CANVAS_N = int(os.environ.get("CANVAS_N", "2400"))
BAR_W = 48
PAD = 45
LOG_VREF = 0.025
PALETTE = np.rint(colormaps["coolwarm"](np.linspace(0, 1, 32))[:, :3] * 255).astype("uint8")


def velocity_fraction(speed: np.ndarray, vmax: float) -> np.ndarray:
    return np.log1p(np.maximum(speed, 0) / LOG_VREF) / np.log1p(vmax / LOG_VREF)


def color_velocity(speed: np.ndarray, vmax: float) -> np.ndarray:
    index = np.minimum((np.clip(velocity_fraction(speed, vmax), 0, 1) * 32).astype(int), 31)
    return PALETTE[index]


def color_maximum(speed: np.ndarray, valid: np.ndarray, field: np.ndarray,
                  mask: np.ndarray, override: str | None = None) -> float:
    if not valid.any() or not np.isfinite(speed[valid]).all():
        raise ValueError("Invalid velocity samples")
    maximum = max(float(speed[valid].max()), float(field[mask > .5].max()), 1e-6)
    if override is not None:
        shared = float(override)
        if not math.isfinite(shared) or shared < maximum:
            raise ValueError("VMAX would hide velocity values; use the full slice maximum")
        maximum = shared
    return maximum


def polygons(geometry: dict) -> list:
    if geometry["type"] == "Polygon":
        return [geometry["coordinates"]]
    if geometry["type"] == "MultiPolygon":
        return geometry["coordinates"]
    raise ValueError(f"Expected polygon geometry, got {geometry['type']}")


def building_mask(features: list, center: tuple, radius: float) -> Image.Image:
    """Union footprints while preserving each polygon's courtyard holes."""
    mask = Image.new("L", (CANVAS_N, CANVAS_N), 0)
    for feature in features:
        for polygon in polygons(feature["geometry"]):
            rings = [[((x - center[0] + radius) / (2 * radius) * (CANVAS_N - 1),
                       (radius - y + center[1]) / (2 * radius) * (CANVAS_N - 1))
                      for x, y, *_ in ring] for ring in polygon]
            xy = np.asarray(rings[0])
            x0, y0 = np.maximum(np.floor(xy.min(axis=0)).astype(int) - 1, 0)
            x1, y1 = np.minimum(np.ceil(xy.max(axis=0)).astype(int) + 2, CANVAS_N)
            if x1 <= x0 or y1 <= y0:
                continue
            tile = Image.new("L", (int(x1 - x0), int(y1 - y0)), 0)
            draw = ImageDraw.Draw(tile)
            for i, ring in enumerate(rings):
                draw.polygon([(x - x0, y - y0) for x, y in ring], fill=255 if i == 0 else 0)
            box = (int(x0), int(y0), int(x1), int(y1))
            mask.paste(ImageChops.lighter(mask.crop(box), tile), box)
    return mask


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
    # Bicubic interpolation can undershoot; a speed cannot be negative.
    return np.maximum(smooth, 0), mask_big


def add_colorbar(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    x0: int,
    y0: int,
    h: int,
    vmax: float,
) -> None:
    vals = LOG_VREF * np.expm1(np.linspace(np.log1p(vmax / LOG_VREF), 0, h))
    bar = color_velocity(vals, vmax).reshape(h, 1, 3)
    bar = np.repeat(bar, BAR_W, axis=1)
    canvas.paste(Image.fromarray(bar, "RGB"), (x0, y0))
    draw.rectangle((x0, y0, x0 + BAR_W, y0 + h), outline=(20, 20, 20), width=1)
    font = ImageFont.truetype(font_manager.findfont("Liberation Sans"), 30)
    last_y = -100
    ticks = sorted(set([0, .02, .05, .1, .2, .5, 1, 2, 4, vmax]), reverse=True)
    for value in ticks:
        if value > vmax:
            continue
        y = y0 + int((1 - velocity_fraction(value, vmax)) * (h - 1))
        if y - last_y < 32:
            continue
        draw.line((x0 + BAR_W, y, x0 + BAR_W + 10, y), fill=(30, 30, 30), width=2)
        draw.text((x0 + BAR_W + 20, y - 14), f"{value:.3g}", font=font, fill=(30, 30, 30))
        last_y = y
    draw.text((x0, y0 - 50), "m/s", font=font, fill=(30, 30, 30))


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
    # The successful mesh is 180 degrees; its original metadata incorrectly said 220.
    arc_width = float(os.environ.get("INLET_WIDTH_DEG", "180"))

    sample = np.load(sample_path)
    speed = sample["speed"].astype(np.float32)
    disk = sample["disk"].astype(bool)
    valid = sample["valid"].astype(bool)
    field, mask = upscale_field(speed, disk, valid)

    vmax = color_maximum(speed, disk & valid, field, mask, os.environ.get("VMAX"))
    rgb = color_velocity(field, vmax)
    alpha = np.clip(mask * 255.0, 0, 255).astype(np.uint8)

    water = np.array([238, 242, 244], dtype=np.uint8)
    blended = (rgb.astype(np.float32) * (alpha[..., None] / 255.0) + water * (1.0 - alpha[..., None] / 255.0)).astype(np.uint8)

    canvas_w = CANVAS_N + 330
    canvas_h = CANVAS_N + 2 * PAD
    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
    canvas.paste(Image.fromarray(blended, "RGB"), (PAD, PAD))
    draw = ImageDraw.Draw(canvas, "RGBA")

    land_path = generated / "sodermalm_osm_island_epsg3006.geojson"
    if land_path.exists():
        land = load_geojson(land_path)["features"][0]["geometry"]
        for ring in polygon_boundary_rings(land):
            draw.line(global_ring_to_px(ring, center, radius), fill=(70, 80, 85), width=2, joint="curve")

    buildings_path = generated / "sodermalm_buildings.geojson"
    if buildings_path.exists() and os.environ.get("SHOW_BUILDINGS", "1") == "1":
        footprints = building_mask(load_geojson(buildings_path)["features"], center, radius)
        canvas.paste((15, 18, 22), (PAD, PAD, PAD + CANVAS_N, PAD + CANVAS_N), footprints)

    if os.environ.get("SHOW_INLET", "1") == "1":
        draw_inlet_arc(draw, radius, wind_from, arc_width)

    if os.environ.get("SHOW_COLORBAR", "1") == "1":
        add_colorbar(canvas, draw, CANVAS_N + 2 * PAD, PAD + CANVAS_N // 8, 3 * CANVAS_N // 4, vmax)

    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out)
    out.with_suffix(".json").write_text(json.dumps({
        "time": float(sample["time"]), "vmax": vmax, "log_reference": LOG_VREF,
        "colormap": "Cool to Warm", "color_bands": 32,
        "inlet_width_degrees": arc_width,
        "buildings": "Opaque footprints; courtyards preserved",
        "sampling": "Nearest-height GLL samples with bicubic display interpolation, not spectral interpolation",
    }, indent=2) + "\n")
    print(out)


if __name__ == "__main__":
    main()
