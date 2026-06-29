#!/usr/bin/env python3
from __future__ import annotations

import os
import struct
import subprocess
from argparse import ArgumentParser
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import PowerNorm
from pymech.neksuite.field import read_header
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt, gaussian_filter

import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import render_urban_depth_scene as scene  # noqa: E402

RUN = Path()
GEOM_FIELD = Path()
FIELD_FILES: list[Path] = []
FRAME_DIR = Path()
OUT = Path()

GRID_N = 320
CHUNK_ELEMS = 160_000
LIFT = 10.0
SMOOTH_SIGMA = 0.75
VMAX_PERCENTILE = 94.0
POWER_GAMMA = 0.45


def field_layout(path: Path):
    with path.open("rb") as f:
        h = read_header(f)
        tag = f.read(4)
        if round(struct.unpack("<f", tag)[0], 5) == 6.54321:
            endian = "<"
        elif round(struct.unpack(">f", tag)[0], 5) == 6.54321:
            endian = ">"
        else:
            raise RuntimeError(f"Could not determine Nek field endianness for {path}")
        f.seek(4 * h.nb_elems_file, os.SEEK_CUR)
        data_start = f.tell()
    bytes_elem = int(h.nb_pts_elem) * h.wdsz
    geom_start = data_start if h.nb_vars[0] else None
    vel_start = data_start + h.nb_elems_file * h.nb_vars[0] * bytes_elem
    return h, endian, bytes_elem, geom_start, vel_start


def fill_missing(grid: np.ndarray) -> np.ndarray:
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


def sample_terrain_following(velocity_field: Path, terrain_interp, xlim, ylim):
    gh, gendian, gbytes_elem, geom_start, _ = field_layout(GEOM_FIELD)
    vh, vendian, vbytes_elem, _, vel_start = field_layout(velocity_field)
    if gh.nb_elems_file != vh.nb_elems_file or gh.nb_pts_elem != vh.nb_pts_elem:
        raise RuntimeError(f"{velocity_field} does not match geometry field layout")

    gdtype = np.dtype(gendian + gh.realtype)
    vdtype = np.dtype(vendian + vh.realtype)
    pts = int(gh.nb_pts_elem)
    elems = int(gh.nb_elems_file)
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

    with GEOM_FIELD.open("rb") as gf, velocity_field.open("rb") as vf:
        for start in range(0, elems, CHUNK_ELEMS):
            stop = min(start + CHUNK_ELEMS, elems)
            ne = stop - start
            geom_count = ne * gh.nb_vars[0] * pts
            vel_count = ne * vh.nb_vars[1] * pts
            gf.seek(geom_start + start * gh.nb_vars[0] * gbytes_elem)
            geom = np.fromfile(gf, dtype=gdtype, count=geom_count).reshape(ne, gh.nb_vars[0], pts)
            vf.seek(vel_start + start * vh.nb_vars[1] * vbytes_elem)
            vel = np.fromfile(vf, dtype=vdtype, count=vel_count).reshape(ne, vh.nb_vars[1], pts)

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
    faces = scene.read_binary_stl(scene.BUILDINGS_FILE, max_faces=200_000)
    normals = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0])
    lengths = np.linalg.norm(normals, axis=1)
    good = lengths > 1.0e-12
    normals[good] /= lengths[good, None]
    roofs = good & (normals[:, 2] > 0.65)
    return [tri[:, :2] for tri in faces[roofs]]


def field_time(path: Path) -> float:
    with path.open("rb") as f:
        return float(read_header(f).time)


def parse_args():
    parser = ArgumentParser(
        description="Render a top-down velocity-magnitude movie at terrain + lift from an immersed_urban_stl run."
    )
    parser.add_argument("run_dir", type=Path, help="Run directory containing fields/ and urban STL symlinks/files.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output mp4 path. Defaults to RUN_DIR/topdown_velocity_magnitude_enhanced.mp4.",
    )
    parser.add_argument(
        "--frame-dir",
        type=Path,
        default=None,
        help="Directory for rendered PNG frames. Defaults to RUN_DIR/movie_frames_enhanced.",
    )
    parser.add_argument("--framerate", type=str, default="2", help="ffmpeg frame rate.")
    parser.add_argument("--vmax-percentile", type=float, default=VMAX_PERCENTILE)
    parser.add_argument("--gamma", type=float, default=POWER_GAMMA)
    return parser.parse_args()


def configure_paths(args) -> None:
    global RUN, GEOM_FIELD, FIELD_FILES, FRAME_DIR, OUT
    RUN = args.run_dir.resolve()
    GEOM_FIELD = RUN / "fields" / "field0.f00000"
    FIELD_FILES = sorted(p for p in (RUN / "fields").glob("field0.f*") if p.suffix != ".nek5000")
    FRAME_DIR = (args.frame_dir or (RUN / "movie_frames_enhanced")).resolve()
    OUT = (args.output or (RUN / "topdown_velocity_magnitude_enhanced.mp4")).resolve()
    scene.TERRAIN_FILE = RUN / "urban_terrain_surface.stl"
    scene.BUILDINGS_FILE = RUN / "urban_buildings.stl"


def main() -> None:
    args = parse_args()
    configure_paths(args)
    FRAME_DIR.mkdir(exist_ok=True)
    terrain_faces = scene.read_binary_stl(scene.TERRAIN_FILE, max_faces=200_000)
    terrain_x, terrain_y, terrain_z = scene.terrain_grid_from_stl(terrain_faces)
    terrain_interp = RegularGridInterpolator(
        (terrain_y[:, 0], terrain_x[0, :]),
        terrain_z,
        bounds_error=False,
        fill_value=float(np.nanmin(terrain_z)),
    )
    xlim = (float(np.nanmin(terrain_x)), float(np.nanmax(terrain_x)))
    ylim = (float(np.nanmin(terrain_y)), float(np.nanmax(terrain_y)))
    roofs = building_roof_polys()

    sampled = []
    vmax_values = []
    for field in FIELD_FILES:
        print(f"Sampling {field.name}", flush=True)
        x, y, u, v = sample_terrain_following(field, terrain_interp, xlim, ylim)
        speed = np.sqrt(u * u + v * v)
        sampled.append((field, x, y, u, v, speed))
        vmax_values.append(float(np.nanpercentile(speed, args.vmax_percentile)))
    vmax = max(vmax_values) if vmax_values else 1.0

    for index, (field, x, y, u, v, speed) in enumerate(sampled):
        t = field_time(field)
        frame = FRAME_DIR / f"frame_{index:04d}.png"
        print(f"Rendering {frame.name} at t={t:g}", flush=True)
        fig, ax = plt.subplots(figsize=(10.5, 10.5), constrained_layout=True)
        im = ax.imshow(
            speed,
            origin="lower",
            extent=[x.min(), x.max(), y.min(), y.max()],
            cmap="turbo",
            norm=PowerNorm(gamma=args.gamma, vmin=0.0, vmax=vmax),
        )
        ax.streamplot(x, y, u, v, color="white", density=2.0, linewidth=0.55, arrowsize=0.45, minlength=0.12)
        ax.streamplot(x, y, u, v, color="black", density=2.0, linewidth=0.18, arrowsize=0.0, minlength=0.12)
        if roofs:
            ax.add_collection(
                PolyCollection(roofs, facecolor="#b7b7b2", edgecolor="#707070", linewidth=0.08, alpha=1.0, zorder=5)
            )
        ax.set_title(f"Enhanced velocity magnitude at terrain + 10 m, t = {t:.3f}")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="|u|")
        fig.savefig(frame, dpi=180)
        plt.close(fig)

    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-framerate",
            args.framerate,
            "-i",
            str(FRAME_DIR / "frame_%04d.png"),
            "-vf",
            "pad=ceil(iw/2)*2:ceil(ih/2)*2",
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            str(OUT),
        ],
        check=True,
    )
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
