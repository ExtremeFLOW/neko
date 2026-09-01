#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import os
import struct
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CASE_ROOT = HERE.parent
GENERATED = Path(os.environ.get("GENERATED_PATH", CASE_ROOT / "generated_2x_xy"))
FIELD = Path(os.environ.get("FIELD_PATH", HERE / "fields" / "field0.f00000"))
GEOMETRY_FIELD = Path(os.environ.get("GEOMETRY_FIELD", HERE / "fields" / "field0.f00000"))
MASK_FIELD = os.environ.get("MASK_FIELD_PATH")
MASK_THRESHOLD = float(os.environ.get("MASK_THRESHOLD", "0.1"))
OUT = Path(os.environ.get("OUT_PATH", HERE / "renders" / "sodermalm_velocity_terrain_plus10.png"))
NPZ_OUT = os.environ.get("NPZ_OUT")
SAMPLE_ONLY = os.environ.get("SAMPLE_ONLY", "0") == "1"
SHOW_STREAMLINES = os.environ.get("SHOW_STREAMLINES", "1") != "0"
VMAX_OVERRIDE = os.environ.get("VMAX")
VMIN_OVERRIDE = os.environ.get("VMIN")
COLOR_SCALE = os.environ.get("COLOR_SCALE", "linear")
COLORMAP_STYLE = os.environ.get("COLORMAP_STYLE", "urban_flow")
LOG_VREF = float(os.environ.get("LOG_VREF", "0.08"))
BUILDING_ALPHA = int(os.environ.get("BUILDING_ALPHA", "185"))
GRID_N = 420
LIFT_M = float(os.environ.get("LIFT_M", "10.0"))
SAMPLE_Z = os.environ.get("SAMPLE_Z")
SAMPLE_TERRAIN_FROM_GEOJSON = os.environ.get("SAMPLE_TERRAIN_FROM_GEOJSON", "0") == "1"
CHUNK_ELEMS = 16000
STREAMLINE_MIN_SPEED = 0.06
STREAMLINE_STEP_M = 24.0
STREAMLINE_MAX_STEPS = 260
STREAMLINE_SEED_SPACING = 34

sys.path.insert(0, str(CASE_ROOT))
from build_sodermalm_cylinder import (  # noqa: E402
    load_geojson,
    polygon_boundary_rings,
    polygon_rings,
    terrain_height,
    terrain_samples,
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


def update_best_scalar(flat, score, value, best_score, best_value) -> None:
    keep = np.isfinite(score) & np.isfinite(value)
    if not np.any(keep):
        return
    flat = flat[keep]
    score = score[keep]
    value = value[keep]
    local = np.full_like(best_score, np.inf)
    np.minimum.at(local, flat, score)
    winners = score <= local[flat]
    flat = flat[winners]
    score = score[winners]
    value = value[winners]
    replace = score < best_score[flat]
    flat = flat[replace]
    best_score[flat] = score[replace]
    best_value[flat] = value[replace]


def read_vector_chunk(handle, offset, start, ne, nelv, npts, dtype, bytes_elem) -> np.ndarray:
    handle.seek(offset + start * 3 * bytes_elem)
    return np.fromfile(handle, dtype=dtype, count=ne * 3 * npts).reshape(ne, 3, npts)


def read_scalar_chunk(handle, offset, scalar_index, start, ne, nelv, npts, dtype, bytes_elem) -> np.ndarray:
    handle.seek(offset + scalar_index * nelv * bytes_elem + start * bytes_elem)
    return np.fromfile(handle, dtype=dtype, count=ne * npts)


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


def velocity_fraction(speed: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    shifted = np.clip(speed - vmin, 0.0, None)
    span = max(vmax - vmin, 1.0e-12)
    if COLOR_SCALE == "log":
        return np.log1p(shifted / LOG_VREF) / np.log1p(span / LOG_VREF)
    return shifted / span


def color_velocity(speed: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    if COLORMAP_STYLE == "urban_flow":
        # Blue-to-red ramp tuned for weak urban-flow variations near 0-6 m/s.
        stops = np.array(
            [
                [54, 44, 133],
                [39, 91, 190],
                [39, 169, 212],
                [53, 207, 161],
                [142, 216, 101],
                [229, 218, 77],
                [241, 152, 55],
                [194, 51, 45],
            ],
            dtype=np.float32,
        )
    else:
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
    t = np.clip(velocity_fraction(speed, vmin, vmax), 0.0, 1.0)
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


def bilinear(grid: np.ndarray, px: float, py: float) -> float:
    if px < 0.0 or py < 0.0 or px >= grid.shape[1] - 1 or py >= grid.shape[0] - 1:
        return float("nan")
    ix = int(math.floor(px))
    iy = int(math.floor(py))
    tx = px - ix
    ty = py - iy
    return float(
        (1.0 - tx) * (1.0 - ty) * grid[iy, ix]
        + tx * (1.0 - ty) * grid[iy, ix + 1]
        + (1.0 - tx) * ty * grid[iy + 1, ix]
        + tx * ty * grid[iy + 1, ix + 1]
    )


def pixel_to_local(px: float, py: float, radius: float) -> tuple[float, float]:
    x = -radius + (px + 0.5) * (2.0 * radius / GRID_N)
    y = radius - (py + 0.5) * (2.0 * radius / GRID_N)
    return x, y


def local_to_pixel(x: float, y: float, radius: float) -> tuple[float, float]:
    px = (x + radius) / (2.0 * radius) * GRID_N - 0.5
    py = (radius - y) / (2.0 * radius) * GRID_N - 0.5
    return px, py


def integrate_streamline(
    seed_px: float,
    seed_py: float,
    u: np.ndarray,
    v: np.ndarray,
    speed: np.ndarray,
    disk: np.ndarray,
    radius: float,
    direction: float,
) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    x, y = pixel_to_local(seed_px, seed_py, radius)
    for _ in range(STREAMLINE_MAX_STEPS):
        px, py = local_to_pixel(x, y, radius)
        ix = int(round(px))
        iy = int(round(py))
        if ix < 0 or iy < 0 or ix >= GRID_N or iy >= GRID_N or not disk[iy, ix]:
            break
        local_speed = bilinear(speed, px, py)
        if not np.isfinite(local_speed) or local_speed < STREAMLINE_MIN_SPEED:
            break
        ux = bilinear(u, px, py)
        vy = bilinear(v, px, py)
        mag = math.hypot(ux, vy)
        if not np.isfinite(mag) or mag < STREAMLINE_MIN_SPEED:
            break
        points.append((x, y))
        x += direction * STREAMLINE_STEP_M * ux / mag
        y += direction * STREAMLINE_STEP_M * vy / mag
    return points


def draw_streamlines(
    draw: ImageDraw.ImageDraw,
    u: np.ndarray,
    v: np.ndarray,
    speed: np.ndarray,
    disk: np.ndarray,
    center: tuple[float, float],
    radius: float,
) -> None:
    for py in range(STREAMLINE_SEED_SPACING // 2, GRID_N, STREAMLINE_SEED_SPACING):
        for px in range(STREAMLINE_SEED_SPACING // 2, GRID_N, STREAMLINE_SEED_SPACING):
            if not disk[py, px] or speed[py, px] < STREAMLINE_MIN_SPEED:
                continue
            backward = integrate_streamline(px, py, u, v, speed, disk, radius, -1.0)
            forward = integrate_streamline(px, py, u, v, speed, disk, radius, 1.0)
            path = list(reversed(backward[1:])) + forward
            if len(path) < 5:
                continue
            global_path = [(x + center[0], y + center[1]) for x, y in path]
            draw_polyline(draw, global_path, center, radius, GRID_N, (245, 248, 250, 150), width=1)

            end = path[-1]
            prev = path[-4] if len(path) >= 4 else path[0]
            dx = end[0] - prev[0]
            dy = end[1] - prev[1]
            norm = math.hypot(dx, dy)
            if norm <= 1.0:
                continue
            ux, uy = dx / norm, dy / norm
            left = (-uy, ux)
            arrow = [
                (end[0] + center[0], end[1] + center[1]),
                (end[0] + center[0] - 22.0 * ux + 8.0 * left[0], end[1] + center[1] - 22.0 * uy + 8.0 * left[1]),
                (end[0] + center[0] - 22.0 * ux - 8.0 * left[0], end[1] + center[1] - 22.0 * uy - 8.0 * left[1]),
            ]
            rows = [project_global(x, y, center, radius, GRID_N) for x, y in arrow]
            draw.polygon(rows, fill=(245, 248, 250, 135))


def inlet_arc_points(
    center: tuple[float, float],
    radius: float,
    wind_from_deg: float,
    arc_width_deg: float,
    n: int = 240,
) -> list[tuple[float, float]]:
    start = wind_from_deg - 0.5 * arc_width_deg
    stop = wind_from_deg + 0.5 * arc_width_deg
    points = []
    for bearing in np.linspace(start, stop, n):
        theta = math.radians(bearing)
        points.append((center[0] + radius * math.sin(theta), center[1] + radius * math.cos(theta)))
    return points


def render_overlay(
    rgb: np.ndarray,
    speed: np.ndarray,
    valid: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    disk: np.ndarray,
    center,
    radius,
    shoreline,
    buildings,
    wind_from_deg: float,
    arc_width_deg: float,
    sample_label: str,
    vmin: float,
    vmax: float,
    time: float,
) -> Image.Image:
    from PIL import Image, ImageDraw

    img = Image.fromarray(rgb, "RGB")
    yy, xx = np.mgrid[0:GRID_N, 0:GRID_N]
    x = -radius + (xx + 0.5) * (2 * radius / GRID_N)
    y = radius - (yy + 0.5) * (2 * radius / GRID_N)
    outside = x * x + y * y > radius * radius
    rgb[outside] = np.array([166, 213, 232], dtype=np.uint8)
    rgb[(~valid) & (~outside)] = np.array([226, 228, 225], dtype=np.uint8)
    img = Image.fromarray(rgb, "RGB")
    draw = ImageDraw.Draw(img, "RGBA")
    if SHOW_STREAMLINES:
        draw_streamlines(draw, u, v, speed, disk & valid, center, radius)

    for ring in shoreline:
        draw_polyline(draw, ring, center, radius, GRID_N, (20, 29, 34, 215), width=2)
    for ring in buildings:
        pts = ring[:-1] if ring and ring[0] == ring[-1] else ring
        if len(pts) < 3:
            continue
        rows = [project_global(x, y, center, radius, GRID_N) for x, y in pts]
        draw.polygon(rows, fill=(185, 185, 178, BUILDING_ALPHA), outline=(82, 82, 76, 170))

    draw_polyline(
        draw,
        inlet_arc_points(center, radius, wind_from_deg, arc_width_deg),
        center,
        radius,
        GRID_N,
        (210, 22, 32, 255),
        width=5,
    )

    canvas = Image.new("RGB", (GRID_N + 92, GRID_N + 36), "white")
    canvas.paste(img, (18, 18))
    cdraw = ImageDraw.Draw(canvas, "RGBA")
    bar_x = GRID_N + 36
    bar_y = 50
    bar_h = GRID_N - 92
    bar_w = 20
    values = np.linspace(vmax, vmin, bar_h, dtype=np.float32)[:, None]
    bar = color_velocity(values, vmin, vmax).reshape(bar_h, 1, 3)
    bar = np.repeat(bar, bar_w, axis=1)
    canvas.paste(Image.fromarray(bar, "RGB"), (bar_x, bar_y))
    cdraw.rectangle((bar_x, bar_y, bar_x + bar_w, bar_y + bar_h), outline=(20, 20, 20, 255), width=1)
    vmid = 0.5 * (vmin + vmax)
    for frac, label in ((0.0, f"{vmax:.2f}"), (0.5, f"{vmid:.2f}"), (1.0, f"{vmin:.2f}")):
        y = bar_y + int(round(frac * bar_h))
        cdraw.line((bar_x + bar_w, y, bar_x + bar_w + 5, y), fill=(20, 20, 20, 255), width=1)
        cdraw.text((bar_x + bar_w + 8, y - 6), label, fill=(20, 20, 20, 255))
    scale_note = "log scale" if COLOR_SCALE == "log" else "linear scale"
    note = f"|u| at {sample_label}, t = {time:.3f}; {scale_note}; inlet {wind_from_deg:.0f} deg +/- {0.5 * arc_width_deg:.0f} deg"
    cdraw.text((18, GRID_N + 20), note, fill=(20, 20, 20, 255))
    cdraw.text((bar_x - 5, bar_y - 22), "|u|", fill=(20, 20, 20, 255))
    return canvas


def main() -> None:
    meta_path = GENERATED / "sodermalm_cylinder_2x_xy_metadata.json"
    if not meta_path.exists():
        meta_path = GENERATED / "sodermalm_cylinder_shallow_debug_metadata.json"
    if not meta_path.exists():
        meta_path = GENERATED / "sodermalm_cylinder_metadata.json"
    if not meta_path.exists():
        meta_path = GENERATED / "sodermalm_cylinder_lumi_coarse_metadata.json"
    meta = json.loads(meta_path.read_text())
    center = tuple(float(v) for v in meta["center_epsg3006"])
    radius = float(meta["radius_m"])
    wind_from_deg = float(meta.get("inflow_from_degrees", 315.0))
    arc_width_deg = float(meta.get("inflow_arc_width_degrees", 90.0))
    land_geometry = load_geojson(GENERATED / "sodermalm_osm_island_epsg3006.geojson")["features"][0]["geometry"]
    shoreline = polygon_boundary_rings(land_geometry)
    terrain_sample_rows = None
    water_level = None
    if SAMPLE_TERRAIN_FROM_GEOJSON:
        terrain_sample_rows, water_level = terrain_samples(GENERATED / "sodermalm_cylinder_contours.geojson", center)

    header = parse_fld_header(FIELD)
    offsets = field_offsets(header)
    if "U" not in offsets:
        raise RuntimeError(f"Field file does not contain U block; rdcode={header['rdcode']!r}")
    geometry_field = FIELD if "X" in offsets else GEOMETRY_FIELD
    geometry_header = parse_fld_header(geometry_field)
    geometry_offsets = field_offsets(geometry_header)
    if "X" not in geometry_offsets:
        raise RuntimeError(f"Geometry field does not contain X block; rdcode={geometry_header['rdcode']!r}")
    for key in ("lx", "ly", "lz", "nelv"):
        if geometry_header[key] != header[key]:
            raise RuntimeError(
                f"Geometry field {geometry_field} is incompatible with {FIELD}: "
                f"{key}={geometry_header[key]} vs {header[key]}"
            )
    dtype = np.dtype(header["endian"] + ("f4" if header["wdsz"] == 4 else "f8"))
    geometry_dtype = np.dtype(geometry_header["endian"] + ("f4" if geometry_header["wdsz"] == 4 else "f8"))
    npts = header["lx"] * header["ly"] * header["lz"]
    nelv = header["nelv"]
    bytes_elem = npts * header["wdsz"]
    geometry_bytes_elem = npts * geometry_header["wdsz"]
    cells = GRID_N * GRID_N
    bottom_z = np.full(cells, np.inf, dtype=np.float32)

    with geometry_field.open("rb") as handle:
        for start in range(0, nelv, CHUNK_ELEMS):
            stop = min(start + CHUNK_ELEMS, nelv)
            ne = stop - start
            xyz = read_vector_chunk(handle, geometry_offsets["X"], start, ne, nelv, npts, geometry_dtype, geometry_bytes_elem)
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
    if SAMPLE_Z is None and SAMPLE_TERRAIN_FROM_GEOJSON:
        terrain = np.full_like(xg, np.nan, dtype=np.float32)
        assert terrain_sample_rows is not None and water_level is not None
        for iy in range(GRID_N):
            for ix in range(GRID_N):
                if disk[iy, ix]:
                    terrain[iy, ix] = terrain_height(
                        float(xg[iy, ix]),
                        float(yg[iy, ix]),
                        terrain_sample_rows,
                        water_level,
                        shoreline,
                        center,
                        radius,
                    )
        target = np.where(disk & (terrain >= bottom - 1.0), terrain + LIFT_M, np.nan).astype(np.float32)
        sample_label = f"GIS terrain + {LIFT_M:g} m"
    elif SAMPLE_Z is None:
        target = np.where(disk, bottom + LIFT_M, np.nan).astype(np.float32)
        sample_label = f"terrain + {LIFT_M:g} m"
    else:
        sample_z = float(SAMPLE_Z)
        target = np.where(disk & (sample_z >= bottom), sample_z, np.nan).astype(np.float32)
        sample_label = f"z = {sample_z:g} m"

    best_score = np.full(cells, np.inf, dtype=np.float32)
    best_u = np.full(cells, np.nan, dtype=np.float32)
    best_v = np.full(cells, np.nan, dtype=np.float32)
    best_w = np.full(cells, np.nan, dtype=np.float32)
    target_flat = target.ravel()

    with geometry_field.open("rb") as geometry_handle, FIELD.open("rb") as field_handle:
        for start in range(0, nelv, CHUNK_ELEMS):
            stop = min(start + CHUNK_ELEMS, nelv)
            ne = stop - start
            xyz = read_vector_chunk(geometry_handle, geometry_offsets["X"], start, ne, nelv, npts, geometry_dtype, geometry_bytes_elem)
            uvw = read_vector_chunk(field_handle, offsets["U"], start, ne, nelv, npts, dtype, bytes_elem)
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
    mask_grid = None
    if MASK_FIELD:
        mask_path = Path(MASK_FIELD)
        mask_header = parse_fld_header(mask_path)
        mask_offsets = field_offsets(mask_header)
        if "S" not in mask_offsets:
            raise RuntimeError(f"Mask field does not contain S block; rdcode={mask_header['rdcode']!r}")
        for key in ("lx", "ly", "lz", "nelv"):
            if mask_header[key] != geometry_header[key]:
                raise RuntimeError(
                    f"Mask field {mask_path} is incompatible with {geometry_field}: "
                    f"{key}={mask_header[key]} vs {geometry_header[key]}"
                )
        mask_dtype = np.dtype(mask_header["endian"] + ("f4" if mask_header["wdsz"] == 4 else "f8"))
        mask_bytes_elem = npts * mask_header["wdsz"]
        best_mask_score = np.full(cells, np.inf, dtype=np.float32)
        best_mask = np.full(cells, np.nan, dtype=np.float32)
        with geometry_field.open("rb") as geometry_handle, mask_path.open("rb") as mask_handle:
            for start in range(0, nelv, CHUNK_ELEMS):
                stop = min(start + CHUNK_ELEMS, nelv)
                ne = stop - start
                xyz = read_vector_chunk(geometry_handle, geometry_offsets["X"], start, ne, nelv, npts, geometry_dtype, geometry_bytes_elem)
                scalar = read_scalar_chunk(mask_handle, mask_offsets["S"], 0, start, ne, nelv, npts, mask_dtype, mask_bytes_elem)
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
                update_best_scalar(flat, score, scalar[inside], best_mask_score, best_mask)
                print(f"sampled mask {stop}/{nelv}", flush=True)
        mask_grid = fill_missing_nearest(best_mask.reshape(GRID_N, GRID_N), np.isfinite(best_mask).reshape(GRID_N, GRID_N))
        valid &= mask_grid < MASK_THRESHOLD
    u = fill_missing_nearest(best_u.reshape(GRID_N, GRID_N), valid)
    v = fill_missing_nearest(best_v.reshape(GRID_N, GRID_N), valid)
    w = fill_missing_nearest(best_w.reshape(GRID_N, GRID_N), valid)
    speed = blur3(np.sqrt(u * u + v * v + w * w), passes=1)
    speed[~valid] = np.nan
    if NPZ_OUT:
        Path(NPZ_OUT).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            NPZ_OUT,
            speed=speed,
            u=u,
            v=v,
            w=w,
            disk=disk,
            valid=valid,
            mask=mask_grid if mask_grid is not None else np.full_like(speed, np.nan, dtype=np.float32),
            x=xg,
            y=yg,
            time=header["time"],
        )
        if SAMPLE_ONLY:
            print(f"Wrote {NPZ_OUT}")
            return
    vmin = float(VMIN_OVERRIDE) if VMIN_OVERRIDE else 0.0
    plot_values = speed[disk & valid]
    if plot_values.size == 0:
        raise RuntimeError("No valid fluid samples after applying the mask threshold")
    vmax = float(VMAX_OVERRIDE) if VMAX_OVERRIDE else float(np.nanpercentile(plot_values, 99.0))
    vmax = max(vmax, vmin + 1.0e-6)
    rgb = color_velocity(np.nan_to_num(speed, nan=vmin), vmin, vmax)

    building_rings = []
    for feature in load_geojson(GENERATED / "sodermalm_buildings.geojson")["features"]:
        for ring in polygon_rings(feature["geometry"]):
            building_rings.append(ring)
    img = render_overlay(
        rgb,
        speed,
        valid,
        u,
        v,
        disk,
        center,
        radius,
        shoreline,
        building_rings,
        wind_from_deg,
        arc_width_deg,
        sample_label,
        vmin,
        vmax,
        header["time"],
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
