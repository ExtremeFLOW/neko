#!/usr/bin/env python3
"""Render the 9x urban result without per-column field sorting."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
from matplotlib.colors import Normalize
from pymech.neksuite import readnek
from scipy.ndimage import distance_transform_edt
from scipy.interpolate import RegularGridInterpolator

import render_urban_depth_scene as scene


HERE = Path(__file__).resolve().parent
FIELD_FILE = HERE / "fine_9x" / "run_21550261" / "fields" / "field0.f00000"
OUTPUT_DIR = HERE / "fine_9x" / "run_21550261" / "renders"
MAIN_OUTPUT = OUTPUT_DIR / "urban_9x_final_streamlines.png"
SLICE_DIR = OUTPUT_DIR / "velocity_slices"

CHUNK = 2_000_000


def read_field(path: Path) -> tuple[np.ndarray, ...]:
    data = readnek(str(path))
    xs, ys, zs, us, vs, ws = [], [], [], [], [], []
    for elem in data.elem:
        xs.append(elem.pos[0].ravel())
        ys.append(elem.pos[1].ravel())
        zs.append(elem.pos[2].ravel())
        us.append(elem.vel[0].ravel())
        vs.append(elem.vel[1].ravel())
        ws.append(elem.vel[2].ravel())
    return tuple(np.concatenate(values) for values in (xs, ys, zs, us, vs, ws))


def update_nearest(
    values: np.ndarray,
    distances: np.ndarray,
    best_values: np.ndarray,
    best_distances: np.ndarray,
) -> None:
    keep = np.isfinite(values) & np.isfinite(distances)
    if not np.any(keep):
        return
    values = values[keep]
    distances = distances[keep]
    order = np.argsort(distances)[::-1]
    values = values[order]
    distances = distances[order]
    update = distances <= best_distances[values]
    best_distances[values[update]] = distances[update]


def ground_velocity_grid(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    terrain_interp: RegularGridInterpolator,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    n: int = 180,
    lift: float = scene.GROUND_LIFT,
) -> tuple[np.ndarray, ...]:
    grid_x = np.linspace(xlim[0], xlim[1], n)
    grid_y = np.linspace(ylim[0], ylim[1], n)
    best_dist = np.full(n * n, np.inf)
    best_u = np.full(n * n, np.nan)
    best_v = np.full(n * n, np.nan)

    for start in range(0, x.size, CHUNK):
        stop = min(start + CHUNK, x.size)
        xc = x[start:stop]
        yc = y[start:stop]
        zc = z[start:stop]
        ix = np.floor((xc - xlim[0]) / (xlim[1] - xlim[0]) * (n - 1)).astype(int)
        iy = np.floor((yc - ylim[0]) / (ylim[1] - ylim[0]) * (n - 1)).astype(int)
        inside = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        if not np.any(inside):
            continue
        flat = iy[inside] * n + ix[inside]
        target = terrain_interp(np.column_stack((yc[inside], xc[inside]))) + lift
        dist = np.abs(zc[inside] - target)
        order = np.argsort(dist)[::-1]
        flat = flat[order]
        dist = dist[order]
        uu = u[start:stop][inside][order]
        vv = v[start:stop][inside][order]
        update = dist <= best_dist[flat]
        best_dist[flat[update]] = dist[update]
        best_u[flat[update]] = uu[update]
        best_v[flat[update]] = vv[update]

    u_grid = best_u.reshape(n, n)
    v_grid = best_v.reshape(n, n)
    for grid in (u_grid, v_grid):
        missing = ~np.isfinite(grid)
        if np.any(missing):
            fill = np.nanmean(grid)
            grid[missing] = 0.0 if not np.isfinite(fill) else fill
    return grid_x, grid_y, u_grid, v_grid


def render_main(x: np.ndarray, y: np.ndarray, z: np.ndarray, u: np.ndarray, v: np.ndarray) -> None:
    terrain_faces = scene.read_binary_stl(scene.TERRAIN_FILE, max_faces=200000)
    building_faces = scene.read_binary_stl(scene.BUILDINGS_FILE, max_faces=200000)
    eye, right, up, forward, x_min, y_min, scale, x_offset, y_offset = scene.make_projection(
        [terrain_faces, building_faces],
        scene.WIDTH,
        scene.HEIGHT,
        scene.ELEVATION,
        scene.AZIMUTH,
        scene.Z_SCALE,
    )

    image = np.zeros((scene.HEIGHT, scene.WIDTH, 3), dtype=np.uint8)
    image[:, :] = np.array([238, 244, 247], dtype=np.uint8)
    depth = np.full((scene.HEIGHT, scene.WIDTH), np.inf)

    terrain_z = np.mean(terrain_faces[:, :, 2], axis=1)
    terrain_norm = Normalize(vmin=float(np.nanmin(terrain_z)), vmax=float(np.nanmax(terrain_z)))
    terrain_base = scene.LAND_CMAP(terrain_norm(terrain_z))[:, :3] * 255.0
    terrain_colors, terrain_valid = scene.shade_faces(terrain_faces, terrain_base)
    terrain_projected = scene.project_points(
        terrain_faces.reshape(-1, 3),
        eye,
        right,
        up,
        forward,
        x_min,
        y_min,
        scale,
        x_offset,
        y_offset,
        scene.HEIGHT,
        scene.Z_SCALE,
    ).reshape((-1, 3, 3))
    scene.draw_faces(image, depth, terrain_faces, terrain_projected, terrain_colors, terrain_valid)

    building_base = np.full((len(building_faces), 3), 178.0)
    building_colors, building_valid = scene.shade_faces(building_faces, building_base)
    building_projected = scene.project_points(
        building_faces.reshape(-1, 3),
        eye,
        right,
        up,
        forward,
        x_min,
        y_min,
        scale,
        x_offset,
        y_offset,
        scene.HEIGHT,
        scene.Z_SCALE,
    ).reshape((-1, 3, 3))
    scene.draw_faces(image, depth, building_faces, building_projected, building_colors, building_valid)

    terrain_x, terrain_y, terrain_grid_z = scene.terrain_grid_from_stl(terrain_faces)
    terrain_interp = RegularGridInterpolator(
        (terrain_y[:, 0], terrain_x[0, :]),
        terrain_grid_z,
        bounds_error=False,
        fill_value=float(np.nanmin(terrain_grid_z)),
    )
    xlim = (float(np.nanmin(terrain_x)), float(np.nanmax(terrain_x)))
    ylim = (float(np.nanmin(terrain_y)), float(np.nanmax(terrain_y)))
    plot_x, plot_y, u_plot, v_plot = ground_velocity_grid(x, y, z, u, v, terrain_interp, xlim, ylim)
    u_interp = RegularGridInterpolator((plot_y, plot_x), u_plot, bounds_error=False, fill_value=0.0)
    v_interp = RegularGridInterpolator((plot_y, plot_x), v_plot, bounds_error=False, fill_value=0.0)

    source_x = xlim[0] + 12.0
    source_y = np.linspace(ylim[0] + 25.0, ylim[1] - 25.0, scene.SOURCE_COUNT)
    seeds = np.column_stack(
        (
            np.full(scene.SOURCE_COUNT, source_x),
            source_y,
            terrain_interp(np.column_stack((source_y, np.full(scene.SOURCE_COUNT, source_x))))
            + scene.GROUND_LIFT,
        )
    )
    paths, path_speeds = scene.trace_streamlines(
        seeds,
        u_interp,
        v_interp,
        terrain_interp,
        xlim,
        ylim,
        scene.GROUND_LIFT,
        scene.STREAMLINE_STEPS,
        scene.STREAMLINE_STEP,
    )

    pil_image = Image.fromarray(image)
    draw = ImageDraw.Draw(pil_image)
    speeds = np.concatenate(path_speeds) if path_speeds else np.array([0.0, 1.0])
    speed_norm = Normalize(vmin=float(np.nanpercentile(speeds, 5)), vmax=float(np.nanpercentile(speeds, 95)))
    for path, path_speed in zip(paths, path_speeds):
        projected = scene.project_points(
            path,
            eye,
            right,
            up,
            forward,
            x_min,
            y_min,
            scale,
            x_offset,
            y_offset,
            scene.HEIGHT,
            scene.Z_SCALE,
        )
        for index in range(len(projected) - 1):
            color = tuple(
                (np.array(scene.FLOW_HEAT_CMAP(speed_norm(path_speed[min(index, len(path_speed) - 1)]))[:3]) * 255)
                .astype(np.uint8)
            )
            draw.line(
                [
                    (float(projected[index, 0]), float(projected[index, 1])),
                    (float(projected[index + 1, 0]), float(projected[index + 1, 1])),
                ],
                fill=color,
                width=scene.LINE_WIDTH,
            )
    pil_image.save(MAIN_OUTPUT)
    print(f"Wrote {MAIN_OUTPUT}")


def nearest_plane(
    a: np.ndarray,
    b: np.ndarray,
    value: np.ndarray,
    distance: np.ndarray,
    alim: tuple[float, float],
    blim: tuple[float, float],
    shape: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ny, nx = shape
    best = np.full(ny * nx, np.inf)
    field = np.full(ny * nx, np.nan)
    for start in range(0, a.size, CHUNK):
        stop = min(start + CHUNK, a.size)
        ac = a[start:stop]
        bc = b[start:stop]
        ix = np.floor((ac - alim[0]) / (alim[1] - alim[0]) * (nx - 1)).astype(int)
        iy = np.floor((bc - blim[0]) / (blim[1] - blim[0]) * (ny - 1)).astype(int)
        inside = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        if not np.any(inside):
            continue
        flat = iy[inside] * nx + ix[inside]
        dist = distance[start:stop][inside]
        vals = value[start:stop][inside]
        order = np.argsort(dist)[::-1]
        flat = flat[order]
        dist = dist[order]
        vals = vals[order]
        update = dist <= best[flat]
        best[flat[update]] = dist[update]
        field[flat[update]] = vals[update]
    aa = np.linspace(alim[0], alim[1], nx)
    bb = np.linspace(blim[0], blim[1], ny)
    return aa, bb, field.reshape(ny, nx)


def plot_slice(path: Path, xx: np.ndarray, yy: np.ndarray, values: np.ndarray, title: str, xlabel: str, ylabel: str) -> None:
    missing = ~np.isfinite(values)
    if np.any(missing) and np.any(~missing):
        nearest = distance_transform_edt(missing, return_distances=False, return_indices=True)
        values = values[tuple(nearest)]

    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    vmax = float(np.nanpercentile(values, 98))
    image = ax.imshow(
        values,
        origin="lower",
        extent=[xx.min(), xx.max(), yy.min(), yy.max()],
        cmap="turbo",
        vmin=0.0,
        vmax=vmax,
        aspect="auto",
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.colorbar(image, ax=ax, label="velocity magnitude")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"Wrote {path}")


def terrain_interpolator() -> RegularGridInterpolator:
    terrain_faces = scene.read_binary_stl(scene.TERRAIN_FILE, max_faces=200000)
    terrain_x, terrain_y, terrain_z = scene.terrain_grid_from_stl(terrain_faces)
    return RegularGridInterpolator(
        (terrain_y[:, 0], terrain_x[0, :]),
        terrain_z,
        bounds_error=False,
        fill_value=float(np.nanmin(terrain_z)),
    )


def render_slices(x: np.ndarray, y: np.ndarray, z: np.ndarray, speed: np.ndarray) -> None:
    SLICE_DIR.mkdir(parents=True, exist_ok=True)
    xlim = (float(np.nanmin(x)), float(np.nanmax(x)))
    ylim = (float(np.nanmin(y)), float(np.nanmax(y)))
    zlim = (float(np.nanmin(z)), float(np.nanmax(z)))
    terrain_interp = terrain_interpolator()

    xx, yy, values = nearest_plane(x, y, speed, np.abs(z - 40.0), xlim, ylim, (420, 420))
    x_grid, y_grid = np.meshgrid(xx, yy)
    below_terrain = 40.0 < terrain_interp(np.column_stack((y_grid.ravel(), x_grid.ravel()))).reshape(values.shape)
    values[below_terrain] = 0.0
    plot_slice(SLICE_DIR / "velocity_slice_z040.png", xx, yy, values, "Velocity magnitude at z = 40 m", "x (m)", "y (m)")

    xx, yy, values = nearest_plane(x, y, speed, np.abs(z - 100.0), xlim, ylim, (420, 420))
    x_grid, y_grid = np.meshgrid(xx, yy)
    below_terrain = 100.0 < terrain_interp(np.column_stack((y_grid.ravel(), x_grid.ravel()))).reshape(values.shape)
    values[below_terrain] = 0.0
    plot_slice(SLICE_DIR / "velocity_slice_z100.png", xx, yy, values, "Velocity magnitude at z = 100 m", "x (m)", "y (m)")

    yy, zz, values = nearest_plane(y, z, speed, np.abs(x - 1000.0), ylim, zlim, (260, 420))
    y_grid, z_grid = np.meshgrid(yy, zz)
    terrain = terrain_interp(np.column_stack((y_grid.ravel(), np.full(y_grid.size, 1000.0)))).reshape(values.shape)
    values[z_grid < terrain] = 0.0
    plot_slice(SLICE_DIR / "velocity_slice_x1000.png", yy, zz, values, "Velocity magnitude at x = 1000 m", "y (m)", "z (m)")

    xx, zz, values = nearest_plane(x, z, speed, np.abs(y - 1000.0), xlim, zlim, (260, 420))
    x_grid, z_grid = np.meshgrid(xx, zz)
    terrain = terrain_interp(np.column_stack((np.full(x_grid.size, 1000.0), x_grid.ravel()))).reshape(values.shape)
    values[z_grid < terrain] = 0.0
    plot_slice(SLICE_DIR / "velocity_slice_y1000.png", xx, zz, values, "Velocity magnitude at y = 1000 m", "x (m)", "z (m)")

    images = [Image.open(path) for path in sorted(SLICE_DIR.glob("velocity_slice_*.png"))]
    width = max(image.width for image in images)
    height = max(image.height for image in images)
    sheet = Image.new("RGB", (2 * width, 2 * height), "white")
    for index, image in enumerate(images):
        sheet.paste(image.convert("RGB"), ((index % 2) * width, (index // 2) * height))
    sheet.save(SLICE_DIR / "velocity_slices_contact_sheet.png")
    print(f"Wrote {SLICE_DIR / 'velocity_slices_contact_sheet.png'}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    x, y, z, u, v, w = read_field(FIELD_FILE)
    speed = np.sqrt(u * u + v * v + w * w)
    render_main(x, y, z, u, v)
    render_slices(x, y, z, speed)


if __name__ == "__main__":
    main()
