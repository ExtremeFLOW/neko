#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import os
import struct
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
CASE_ROOT = HERE.parent
GENERATED = CASE_ROOT / "generated_2x_xy"
FIELD = HERE / "fields" / "field0.f00000"
OUT = HERE / "renders" / "sodermalm_velocity_terrain_plus10.png"
GRID_N = 420
LIFT_M = 10.0
CHUNK_ELEMS = 16000

sys.path.insert(0, str(CASE_ROOT))
from build_sodermalm_cylinder import (  # noqa: E402
    load_geojson,
    polygon_boundary_rings,
    polygon_rings,
)


def parse_fld_header(path: Path) -> dict:
    with path.open("rb") as handle:
        raw = handle.read(132)
        header = raw.decode("ascii", errors="replace")
        marker = handle.read(4)
    parts = header.split()
    if len(parts) < 11 or parts[0] != "#std":
        raise RuntimeError(f"Unexpected field header: {header!r}")
    little = round(struct.unpack("<f", marker)[0], 5)
    big = round(struct.unpack(">f", marker)[0], 5)
    if little == 6.54321:
        endian = "<"
    elif big == 6.54321:
        endian = ">"
    else:
        raise RuntimeError("Could not determine field file endianness")
    return {
        "wdsz": int(parts[1]),
        "lx": int(parts[2]),
        "ly": int(parts[3]),
        "lz": int(parts[4]),
        "nelv": int(parts[5]),
        "time": float(parts[7].replace("D", "E")),
        "step": int(parts[8]),
        "rdcode": "".join(parts[11:]),
        "endian": endian,
    }


def field_offsets(header: dict) -> dict[str, int]:
    npts = header["lx"] * header["ly"] * header["lz"]
    nelv = header["nelv"]
    scalar_bytes = npts * header["wdsz"]
    offset = 132 + 4 + 4 * nelv
    offsets: dict[str, int] = {}
    rdcode = header["rdcode"]
    i = 0
    while i < len(rdcode):
        code = rdcode[i]
        if code == "X":
            offsets["X"] = offset
            offset += nelv * 3 * scalar_bytes
        elif code == "U":
            offsets["U"] = offset
            offset += nelv * 3 * scalar_bytes
        elif code in ("P", "T"):
            offsets[code] = offset
            offset += nelv * scalar_bytes
        elif code == "S":
            n_scalars = int(rdcode[i + 1 : i + 3])
            offsets["S"] = offset
            offset += nelv * n_scalars * scalar_bytes
            i += 2
        i += 1
    return offsets


def update_best(flat, score, u, v, w, best_score, best_u, best_v, best_w) -> None:
    keep = np.isfinite(score) & np.isfinite(u) & np.isfinite(v) & np.isfinite(w)
    if not np.any(keep):
        return
    flat = flat[keep]
    score = score[keep]
    u = u[keep]
    v = v[keep]
    w = w[keep]
    local = np.full_like(best_score, np.inf)
    np.minimum.at(local, flat, score)
    winners = score <= local[flat]
    flat = flat[winners]
    score = score[winners]
    u = u[winners]
    v = v[winners]
    w = w[winners]
    replace = score < best_score[flat]
    flat = flat[replace]
    best_score[flat] = score[replace]
    best_u[flat] = u[replace]
    best_v[flat] = v[replace]
    best_w[flat] = w[replace]


def fill_missing_nearest(grid: np.ndarray, valid_mask: np.ndarray) -> np.ndarray:
    if not np.any(valid_mask):
        return np.zeros_like(grid)
    filled = grid.copy()
    valid = valid_mask.copy()
    filled[~valid] = 0.0
    for _ in range(max(grid.shape)):
        if np.all(valid):
            break
        new_filled = filled.copy()
        new_valid = valid.copy()
        count = np.zeros_like(filled, dtype=np.float32)
        acc = np.zeros_like(filled, dtype=np.float32)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dx == 0 and dy == 0:
                    continue
                shifted_valid = np.roll(np.roll(valid, dy, axis=0), dx, axis=1)
                shifted_values = np.roll(np.roll(filled, dy, axis=0), dx, axis=1)
                if dy == -1:
                    shifted_valid[-1, :] = False
                elif dy == 1:
                    shifted_valid[0, :] = False
                if dx == -1:
                    shifted_valid[:, -1] = False
                elif dx == 1:
                    shifted_valid[:, 0] = False
                take = (~valid) & shifted_valid
                acc[take] += shifted_values[take]
                count[take] += 1.0
        grown = (~valid) & (count > 0.0)
        new_filled[grown] = acc[grown] / count[grown]
        new_valid[grown] = True
        filled, valid = new_filled, new_valid
    return filled


def blur3(grid: np.ndarray, passes: int = 2) -> np.ndarray:
    out = grid.astype(np.float32, copy=True)
    for _ in range(passes):
        pad = np.pad(out, 1, mode="edge")
        out = (
            4.0 * pad[1:-1, 1:-1]
            + 2.0 * (pad[:-2, 1:-1] + pad[2:, 1:-1] + pad[1:-1, :-2] + pad[1:-1, 2:])
            + pad[:-2, :-2]
            + pad[:-2, 2:]
            + pad[2:, :-2]
            + pad[2:, 2:]
        ) / 16.0
    return out


def color_velocity(speed: np.ndarray, vmax: float) -> np.ndarray:
    # A restrained viridis-like ramp, kept readable without saturating the whole map.
    stops = np.array(
        [
            [39, 35, 74],
            [48, 103, 141],
            [53, 157, 139],
            [128, 198, 97],
            [235, 220, 77],
            [246, 146, 54],
            [166, 37, 41],
        ],
        dtype=np.float32,
    )
    t = np.clip(speed / max(vmax, 1.0e-12), 0.0, 1.0)
    pos = t * (len(stops) - 1)
    idx = np.floor(pos).astype(np.int32)
    frac = (pos - idx)[..., None]
    idx1 = np.minimum(idx + 1, len(stops) - 1)
    rgb = stops[idx] * (1.0 - frac) + stops[idx1] * frac
    return np.uint8(np.clip(rgb, 0, 255))


def project_global(x: float, y: float, center: tuple[float, float], radius: float, size: int) -> tuple[int, int]:
    px = int(round((x - (center[0] - radius)) / (2 * radius) * (size - 1)))
    py = int(round(((center[1] + radius) - y) / (2 * radius) * (size - 1)))
    return px, py


def draw_polyline(draw: ImageDraw.ImageDraw, points, center, radius, size, fill, width=1) -> None:
    rows = [project_global(x, y, center, radius, size) for x, y in points]
    if len(rows) >= 2:
        draw.line(rows, fill=fill, width=width, joint="curve")


def render_overlay(rgb: np.ndarray, speed: np.ndarray, valid: np.ndarray, center, radius, shoreline, buildings, vmax: float, time: float) -> Image.Image:
    img = Image.fromarray(rgb, "RGB")
    alpha = Image.new("L", (GRID_N, GRID_N), 255)
    alpha_arr = np.asarray(alpha)
    yy, xx = np.mgrid[0:GRID_N, 0:GRID_N]
    x = -radius + (xx + 0.5) * (2 * radius / GRID_N)
    y = radius - (yy + 0.5) * (2 * radius / GRID_N)
    outside = x * x + y * y > radius * radius
    rgb[outside] = np.array([166, 213, 232], dtype=np.uint8)
    img = Image.fromarray(rgb, "RGB")
    draw = ImageDraw.Draw(img, "RGBA")

    for ring in shoreline:
        draw_polyline(draw, ring, center, radius, GRID_N, (20, 29, 34, 215), width=2)
    for ring in buildings:
        pts = ring[:-1] if ring and ring[0] == ring[-1] else ring
        if len(pts) < 3:
            continue
        rows = [project_global(x, y, center, radius, GRID_N) for x, y in pts]
        draw.polygon(rows, fill=(185, 185, 178, 185), outline=(82, 82, 76, 170))

    # Red inlet arc.
    inlet = [
        (center[0] + radius * math.cos(math.radians(a)), center[1] + radius * math.sin(math.radians(a)))
        for a in np.linspace(90.0, 180.0, 120)
    ]
    draw_polyline(draw, inlet, center, radius, GRID_N, (210, 22, 32, 255), width=5)

    canvas = Image.new("RGB", (GRID_N + 92, GRID_N + 36), "white")
    canvas.paste(img, (18, 18))
    cdraw = ImageDraw.Draw(canvas, "RGBA")
    bar_x = GRID_N + 36
    bar_y = 50
    bar_h = GRID_N - 92
    bar_w = 20
    values = np.linspace(vmax, 0.0, bar_h, dtype=np.float32)[:, None]
    bar = color_velocity(values, vmax).reshape(bar_h, 1, 3)
    bar = np.repeat(bar, bar_w, axis=1)
    canvas.paste(Image.fromarray(bar, "RGB"), (bar_x, bar_y))
    cdraw.rectangle((bar_x, bar_y, bar_x + bar_w, bar_y + bar_h), outline=(20, 20, 20, 255), width=1)
    for frac, label in ((0.0, f"{vmax:.2f}"), (0.5, f"{0.5 * vmax:.2f}"), (1.0, "0")):
        y = bar_y + int(round(frac * bar_h))
        cdraw.line((bar_x + bar_w, y, bar_x + bar_w + 5, y), fill=(20, 20, 20, 255), width=1)
        cdraw.text((bar_x + bar_w + 8, y - 6), label, fill=(20, 20, 20, 255))
    cdraw.text((18, GRID_N + 20), f"|u| at terrain + {LIFT_M:g} m, t = {time:.3f}", fill=(20, 20, 20, 255))
    cdraw.text((bar_x - 5, bar_y - 22), "|u|", fill=(20, 20, 20, 255))
    return canvas


def main() -> None:
    meta_path = GENERATED / "sodermalm_cylinder_2x_xy_metadata.json"
    if not meta_path.exists():
        meta_path = GENERATED / "sodermalm_cylinder_metadata.json"
    meta = json.loads(meta_path.read_text())
    center = tuple(float(v) for v in meta["center_epsg3006"])
    radius = float(meta["radius_m"])
    land_geometry = load_geojson(GENERATED / "sodermalm_osm_island_epsg3006.geojson")["features"][0]["geometry"]
    shoreline = polygon_boundary_rings(land_geometry)

    header = parse_fld_header(FIELD)
    offsets = field_offsets(header)
    if "X" not in offsets or "U" not in offsets:
        raise RuntimeError(f"Field file does not contain X and U blocks; rdcode={header['rdcode']!r}")
    dtype = np.dtype(header["endian"] + ("f4" if header["wdsz"] == 4 else "f8"))
    npts = header["lx"] * header["ly"] * header["lz"]
    nelv = header["nelv"]
    bytes_elem = npts * header["wdsz"]
    cells = GRID_N * GRID_N
    bottom_z = np.full(cells, np.inf, dtype=np.float32)

    with FIELD.open("rb") as handle:
        for start in range(0, nelv, CHUNK_ELEMS):
            stop = min(start + CHUNK_ELEMS, nelv)
            ne = stop - start
            handle.seek(offsets["X"] + start * 3 * bytes_elem)
            xyz = np.fromfile(handle, dtype=dtype, count=ne * 3 * npts).reshape(ne, 3, npts)
            x = xyz[:, 0, :].ravel()
            y = xyz[:, 1, :].ravel()
            z = xyz[:, 2, :].ravel()
            ix = np.floor((x + radius) / (2 * radius) * GRID_N).astype(np.int32)
            iy = np.floor((radius - y) / (2 * radius) * GRID_N).astype(np.int32)
            inside = (ix >= 0) & (ix < GRID_N) & (iy >= 0) & (iy < GRID_N)
            flat = iy[inside] * GRID_N + ix[inside]
            np.minimum.at(bottom_z, flat, z[inside].astype(np.float32))
            print(f"bottom scan {stop}/{nelv}", flush=True)

    bottom_valid = np.isfinite(bottom_z).reshape(GRID_N, GRID_N)
    bottom = fill_missing_nearest(bottom_z.reshape(GRID_N, GRID_N), bottom_valid)
    yy, xx = np.mgrid[0:GRID_N, 0:GRID_N]
    xg = -radius + (xx + 0.5) * (2 * radius / GRID_N)
    yg = radius - (yy + 0.5) * (2 * radius / GRID_N)
    disk = xg * xg + yg * yg <= radius * radius
    target = np.where(disk, bottom + LIFT_M, np.nan).astype(np.float32)

    best_score = np.full(cells, np.inf, dtype=np.float32)
    best_u = np.full(cells, np.nan, dtype=np.float32)
    best_v = np.full(cells, np.nan, dtype=np.float32)
    best_w = np.full(cells, np.nan, dtype=np.float32)
    target_flat = target.ravel()

    with FIELD.open("rb") as handle:
        for start in range(0, nelv, CHUNK_ELEMS):
            stop = min(start + CHUNK_ELEMS, nelv)
            ne = stop - start
            handle.seek(offsets["X"] + start * 3 * bytes_elem)
            xyz = np.fromfile(handle, dtype=dtype, count=ne * 3 * npts).reshape(ne, 3, npts)
            handle.seek(offsets["U"] + start * 3 * bytes_elem)
            uvw = np.fromfile(handle, dtype=dtype, count=ne * 3 * npts).reshape(ne, 3, npts)
            x = xyz[:, 0, :].ravel()
            y = xyz[:, 1, :].ravel()
            z = xyz[:, 2, :].ravel()
            ix = np.floor((x + radius) / (2 * radius) * GRID_N).astype(np.int32)
            iy = np.floor((radius - y) / (2 * radius) * GRID_N).astype(np.int32)
            inside = (ix >= 0) & (ix < GRID_N) & (iy >= 0) & (iy < GRID_N)
            if not np.any(inside):
                continue
            flat = iy[inside] * GRID_N + ix[inside]
            score = np.abs(z[inside] - target_flat[flat]).astype(np.float32)
            update_best(flat, score, uvw[:, 0, :].ravel()[inside], uvw[:, 1, :].ravel()[inside], uvw[:, 2, :].ravel()[inside], best_score, best_u, best_v, best_w)
            print(f"sampled {stop}/{nelv}", flush=True)

    valid = np.isfinite(best_u).reshape(GRID_N, GRID_N)
    u = fill_missing_nearest(best_u.reshape(GRID_N, GRID_N), valid)
    v = fill_missing_nearest(best_v.reshape(GRID_N, GRID_N), valid)
    w = fill_missing_nearest(best_w.reshape(GRID_N, GRID_N), valid)
    speed = blur3(np.sqrt(u * u + v * v + w * w), passes=1)
    vmax = float(np.nanpercentile(speed[disk], 99.0))
    vmax = max(vmax, 1.0e-6)
    rgb = color_velocity(speed, vmax)

    building_rings = []
    for feature in load_geojson(GENERATED / "sodermalm_buildings.geojson")["features"]:
        for ring in polygon_rings(feature["geometry"]):
            building_rings.append(ring)
    img = render_overlay(rgb, speed, valid, center, radius, shoreline, building_rings, vmax, header["time"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
