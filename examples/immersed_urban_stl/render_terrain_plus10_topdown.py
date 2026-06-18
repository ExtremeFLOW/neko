#!/usr/bin/env python3
from __future__ import annotations

import os
import struct
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/cfs/klemming/scratch/t/tadao/cpu-smoke-current/cache/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from pymech.neksuite.field import read_header
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt, gaussian_filter

import renders_scene as scene

ROOT = Path("/cfs/klemming/scratch/t/tadao/cpu-smoke-current")
RUN = ROOT / "runs" / "urban-18x-main-cache-21584753"
GEOM = ROOT / "immersed-urban-input-18x" / "generated"
FIELD = RUN / "fields" / "field0.f00000"
OUT = RUN / "renders" / "terrain_plus_10m_topdown.png"

scene.TERRAIN_FILE = GEOM / "terrain.stl"
scene.BUILDINGS_FILE = GEOM / "buildings.stl"

GRID_N = 320
CHUNK_ELEMS = 160_000
LIFT = 10.0
SMOOTH_SIGMA = 0.75


def field_layout(path: Path):
    with path.open("rb") as f:
        h = read_header(f)
        tag = f.read(4)
        if round(struct.unpack("<f", tag)[0], 5) == 6.54321:
            endian = "<"
        elif round(struct.unpack(">f", tag)[0], 5) == 6.54321:
            endian = ">"
        else:
            raise RuntimeError("Could not determine Nek field endianness")
        f.seek(4 * h.nb_elems_file, os.SEEK_CUR)
        data_start = f.tell()
    bytes_elem = int(h.nb_pts_elem) * h.wdsz
    geom_start = data_start
    vel_start = geom_start + h.nb_elems_file * h.nb_vars[0] * bytes_elem
    return h, endian, bytes_elem, geom_start, vel_start


def fill_missing(grid):
    missing = ~np.isfinite(grid)
    if not np.any(missing):
        return grid
    if np.all(missing):
        grid[:] = 0.0
        return grid
    idx = distance_transform_edt(missing, return_distances=False, return_indices=True)
    grid[missing] = grid[tuple(idx[:, missing])]
    return grid


def update_min(flat, score, uu, vv, best_score, best_u, best_v):
    keep = np.isfinite(score) & np.isfinite(uu) & np.isfinite(vv)
    if not np.any(keep):
        return
    flat = flat[keep]
    score = score[keep]
    uu = uu[keep]
    vv = vv[keep]
    local = np.full_like(best_score, np.inf)
    np.minimum.at(local, flat, score)
    winners = score <= local[flat]
    flat = flat[winners]
    score = score[winners]
    uu = uu[winners]
    vv = vv[winners]
    update = score < best_score[flat]
    flat = flat[update]
    best_score[flat] = score[update]
    best_u[flat] = uu[update]
    best_v[flat] = vv[update]


def sample_terrain_following(field, terrain_interp, xlim, ylim):
    h, endian, bytes_elem, geom_start, vel_start = field_layout(field)
    dtype = np.dtype(endian + h.realtype)
    pts = int(h.nb_pts_elem)
    elems = int(h.nb_elems_file)
    cells = GRID_N * GRID_N

    below_d = np.full(cells, np.inf, dtype=np.float32)
    above_d = np.full(cells, np.inf, dtype=np.float32)
    near_d = np.full(cells, np.inf, dtype=np.float32)
    below_u = np.full(cells, np.nan, dtype=np.float32)
    below_v = np.full(cells, np.nan, dtype=np.float32)
    above_u = np.full(cells, np.nan, dtype=np.float32)
    above_v = np.full(cells, np.nan, dtype=np.float32)
    near_u = np.full(cells, np.nan, dtype=np.float32)
    near_v = np.full(cells, np.nan, dtype=np.float32)

    with field.open("rb") as f:
        for start in range(0, elems, CHUNK_ELEMS):
            stop = min(start + CHUNK_ELEMS, elems)
            ne = stop - start
            count = ne * h.nb_vars[0] * pts
            f.seek(geom_start + start * h.nb_vars[0] * bytes_elem)
            geom = np.fromfile(f, dtype=dtype, count=count).reshape(ne, h.nb_vars[0], pts)
            f.seek(vel_start + start * h.nb_vars[1] * bytes_elem)
            vel = np.fromfile(f, dtype=dtype, count=count).reshape(ne, h.nb_vars[1], pts)

            x = geom[:, 0, :].ravel()
            y = geom[:, 1, :].ravel()
            z = geom[:, 2, :].ravel()
            u = vel[:, 0, :].ravel()
            v = vel[:, 1, :].ravel()

            ix = np.floor((x - xlim[0]) / (xlim[1] - xlim[0]) * (GRID_N - 1)).astype(np.int32)
            iy = np.floor((y - ylim[0]) / (ylim[1] - ylim[0]) * (GRID_N - 1)).astype(np.int32)
            inside = (ix >= 0) & (ix < GRID_N) & (iy >= 0) & (iy < GRID_N)
            if not np.any(inside):
                continue

            flat = iy[inside] * GRID_N + ix[inside]
            target = terrain_interp(np.column_stack((y[inside], x[inside]))) + LIFT
            dz = z[inside] - target
            update_min(flat, np.abs(dz), u[inside], v[inside], near_d, near_u, near_v)
            mask = dz <= 0.0
            update_min(flat[mask], -dz[mask], u[inside][mask], v[inside][mask], below_d, below_u, below_v)
            mask = dz >= 0.0
            update_min(flat[mask], dz[mask], u[inside][mask], v[inside][mask], above_d, above_u, above_v)

            if start == 0 or stop == elems or (start // CHUNK_ELEMS) % 10 == 0:
                bracketed = np.count_nonzero(np.isfinite(below_u) & np.isfinite(above_u))
                print(f"sampled {stop}/{elems}; bracketed {bracketed}/{cells}", flush=True)

    bracketed = np.isfinite(below_u) & np.isfinite(above_u)
    weight = below_d / (below_d + above_d)
    u = np.where(bracketed, below_u * (1.0 - weight) + above_u * weight, near_u).reshape(GRID_N, GRID_N)
    v = np.where(bracketed, below_v * (1.0 - weight) + above_v * weight, near_v).reshape(GRID_N, GRID_N)
    u = gaussian_filter(fill_missing(u), SMOOTH_SIGMA, mode="nearest")
    v = gaussian_filter(fill_missing(v), SMOOTH_SIGMA, mode="nearest")
    x = np.linspace(xlim[0], xlim[1], GRID_N)
    y = np.linspace(ylim[0], ylim[1], GRID_N)
    return x, y, u, v


def building_roof_polys():
    faces = scene.read_binary_stl(scene.BUILDINGS_FILE, max_faces=200000)
    normals = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0])
    lengths = np.linalg.norm(normals, axis=1)
    good = lengths > 1.0e-12
    normals[good] /= lengths[good, None]
    roofs = good & (normals[:, 2] > 0.65)
    polys = [tri[:, :2] for tri in faces[roofs]]
    return polys


def main():
    terrain_faces = scene.read_binary_stl(scene.TERRAIN_FILE, max_faces=200000)
    terrain_x, terrain_y, terrain_z = scene.terrain_grid_from_stl(terrain_faces)
    terrain_interp = RegularGridInterpolator(
        (terrain_y[:, 0], terrain_x[0, :]), terrain_z, bounds_error=False, fill_value=float(np.nanmin(terrain_z))
    )
    xlim = (float(np.nanmin(terrain_x)), float(np.nanmax(terrain_x)))
    ylim = (float(np.nanmin(terrain_y)), float(np.nanmax(terrain_y)))
    x, y, u, v = sample_terrain_following(FIELD, terrain_interp, xlim, ylim)
    speed = np.sqrt(u*u + v*v)

    fig, ax = plt.subplots(figsize=(10.5, 10.5), constrained_layout=True)
    norm = Normalize(vmin=0.0, vmax=float(np.nanpercentile(speed, 98.5)))
    im = ax.imshow(speed, origin="lower", extent=[x.min(), x.max(), y.min(), y.max()], cmap="turbo", norm=norm)

    # Streamlines are drawn first, then roofs hide any line segments passing through buildings.
    ax.streamplot(x, y, u, v, color="white", density=2.0, linewidth=0.55, arrowsize=0.45, minlength=0.12)
    ax.streamplot(x, y, u, v, color="black", density=2.0, linewidth=0.18, arrowsize=0.0, minlength=0.12)

    roofs = building_roof_polys()
    if roofs:
        ax.add_collection(PolyCollection(roofs, facecolor="#b7b7b2", edgecolor="#707070", linewidth=0.08, alpha=1.0, zorder=5))

    ax.set_title("Velocity magnitude and streamlines at terrain + 10 m, t = 0.006")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="|u|")
    fig.savefig(OUT, dpi=210)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
