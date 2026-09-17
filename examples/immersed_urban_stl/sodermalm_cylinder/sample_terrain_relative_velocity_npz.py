#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import struct
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CASE_ROOT = HERE
GENERATED = Path(os.environ.get("GENERATED_PATH", HERE / "generated_sharp_mask"))
FIELD = Path(os.environ.get("FIELD_PATH", HERE / "run" / "fields" / "field0.f00000"))
GEOMETRY_FIELD = Path(os.environ.get("GEOMETRY_FIELD", HERE / "run" / "fields" / "field0.f00000"))
NPZ_OUT = Path(os.environ.get("NPZ_OUT", HERE / "renders" / "terrain_plus20.npz"))
GRID_N = 420
LIFT_M = float(os.environ.get("LIFT_M", "20.0"))
SAMPLE_Z = os.environ.get("SAMPLE_Z")
SAMPLE_TERRAIN_FROM_GEOJSON = os.environ.get("SAMPLE_TERRAIN_FROM_GEOJSON", "0") == "1"
CHUNK_ELEMS = 16000
sys.path.insert(0, str(CASE_ROOT))
from build_sodermalm_cylinder import (  # noqa: E402
    load_geojson,
    polygon_boundary_rings,
    terrain_height,
    terrain_samples,
)


def metadata_path(generated: Path) -> Path:
    matches = sorted(p for p in generated.glob("*metadata.json") if not p.name.startswith("._"))
    if matches:
        return matches[0]
    raise FileNotFoundError(f"No metadata JSON found in {generated}")


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


def read_vector_chunk(handle, offset, start, ne, nelv, npts, dtype, bytes_elem) -> np.ndarray:
    handle.seek(offset + start * 3 * bytes_elem)
    return np.fromfile(handle, dtype=dtype, count=ne * 3 * npts).reshape(ne, 3, npts)


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


def main() -> None:
    meta_path = metadata_path(GENERATED)
    meta = json.loads(meta_path.read_text())
    center = tuple(float(v) for v in meta["center_epsg3006"])
    radius = float(meta["radius_m"])
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
    ids = np.fromfile(FIELD, dtype=header["endian"] + "i4", count=nelv, offset=136)
    geometry_ids = np.fromfile(geometry_field, dtype=geometry_header["endian"] + "i4",
                               count=nelv, offset=136)
    if len(ids) != nelv or not np.array_equal(ids, geometry_ids):
        raise ValueError("Velocity and geometry fields have different element ordering")
    if FIELD.stat().st_size < offsets["U"] + nelv * 3 * bytes_elem:
        raise ValueError("Incomplete velocity field")
    if geometry_field.stat().st_size < geometry_offsets["X"] + nelv * 3 * geometry_bytes_elem:
        raise ValueError("Incomplete geometry field")
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
    u = fill_missing_nearest(best_u.reshape(GRID_N, GRID_N), valid)
    v = fill_missing_nearest(best_v.reshape(GRID_N, GRID_N), valid)
    w = fill_missing_nearest(best_w.reshape(GRID_N, GRID_N), valid)
    speed = blur3(np.sqrt(u * u + v * v + w * w), passes=1)
    speed[~valid] = np.nan
    NPZ_OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        NPZ_OUT,
        speed=speed, u=u, v=v, w=w, disk=disk, valid=valid,
        x=xg, y=yg, time=header["time"], sample_label=sample_label,
    )
    print(f"Wrote {NPZ_OUT}")


if __name__ == "__main__":
    main()
