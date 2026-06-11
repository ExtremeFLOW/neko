#!/usr/bin/env python3
"""Depth-buffered full-scene render for the urban immersed-boundary example."""

from __future__ import annotations

import os
import struct
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
from PIL import Image, ImageDraw
from matplotlib.colors import LinearSegmentedColormap, Normalize
from pymech.neksuite import readnek
from scipy.interpolate import RegularGridInterpolator, griddata


LAND_CMAP = LinearSegmentedColormap.from_list(
    "urban_land",
    [
        "#60714a",
        "#7e8758",
        "#999166",
        "#b0a477",
        "#c0b58e",
        "#d0c7aa",
        "#ddd8c7",
    ],
)

FLOW_HEAT_CMAP = LinearSegmentedColormap.from_list(
    "urban_flow_heat",
    [
        "#001f9e",
        "#0056ff",
        "#00b7ff",
        "#20e5c5",
        "#b7f34a",
        "#ffe44a",
        "#ff9c18",
    ],
)

HERE = Path(__file__).resolve().parent
FIELD_FILE = HERE / "fields" / "field0.f00000"
TERRAIN_FILE = HERE / "urban_terrain_surface.stl"
BUILDINGS_FILE = HERE / "urban_buildings.stl"
OUTPUT_FILE = HERE / "urban_depth_full_ground_streamlines_view_restored.png"

WIDTH = 1800
HEIGHT = 1250
SOURCE_COUNT = 56
GROUND_LIFT = 10.0
STREAMLINE_STEP = 8.0
STREAMLINE_STEPS = 520
LINE_WIDTH = 2
ELEVATION = 42.0
AZIMUTH = 318.0
Z_SCALE = 3.0


def read_binary_stl(path: Path, max_faces: int = 10000) -> np.ndarray:
    data = path.read_bytes()
    ntri = struct.unpack_from("<I", data, 80)[0]
    step = max(1, ntri // max_faces)
    faces = []
    offset = 84
    for tri_index in range(ntri):
        offset += 12
        triangle = []
        for _ in range(3):
            triangle.append(struct.unpack_from("<fff", data, offset))
            offset += 12
        offset += 2
        if tri_index % step == 0:
            faces.append(triangle)
    return np.asarray(faces, dtype=float)


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
    return (
        np.concatenate(xs),
        np.concatenate(ys),
        np.concatenate(zs),
        np.concatenate(us),
        np.concatenate(vs),
        np.concatenate(ws),
    )


def terrain_grid_from_stl(faces: np.ndarray) -> tuple[np.ndarray, ...]:
    vertices = faces.reshape(-1, 3)
    unique_xy, inverse = np.unique(vertices[:, :2], axis=0, return_inverse=True)
    terrain_z = np.full(unique_xy.shape[0], np.inf)
    np.minimum.at(terrain_z, inverse, vertices[:, 2])
    nx = min(180, max(40, len(np.unique(unique_xy[:, 0]))))
    ny = min(180, max(40, len(np.unique(unique_xy[:, 1]))))
    xi = np.linspace(float(np.min(unique_xy[:, 0])), float(np.max(unique_xy[:, 0])), nx)
    yi = np.linspace(float(np.min(unique_xy[:, 1])), float(np.max(unique_xy[:, 1])), ny)
    x_grid, y_grid = np.meshgrid(xi, yi)
    z_grid = griddata(unique_xy, terrain_z, (x_grid, y_grid), method="linear")
    if np.isnan(z_grid).any():
        nearest = griddata(unique_xy, terrain_z, (x_grid, y_grid), method="nearest")
        z_grid = np.where(np.isnan(z_grid), nearest, z_grid)
    return x_grid, y_grid, z_grid


def regularize_vector_field(
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    u_grid: np.ndarray,
    v_grid: np.ndarray,
    n: int = 150,
) -> tuple[np.ndarray, ...]:
    xlim = (float(np.nanmin(x_grid)), float(np.nanmax(x_grid)))
    ylim = (float(np.nanmin(y_grid)), float(np.nanmax(y_grid)))
    plot_x = np.linspace(xlim[0], xlim[1], n)
    plot_y = np.linspace(ylim[0], ylim[1], n)
    plot_x_grid, plot_y_grid = np.meshgrid(plot_x, plot_y)
    source = np.column_stack((x_grid.ravel(), y_grid.ravel()))
    u = griddata(source, u_grid.ravel(), (plot_x_grid, plot_y_grid), method="linear")
    v = griddata(source, v_grid.ravel(), (plot_x_grid, plot_y_grid), method="linear")
    for grid, values in ((u, u_grid.ravel()), (v, v_grid.ravel())):
        missing = np.isnan(grid)
        if np.any(missing):
            fill = griddata(source, values, (plot_x_grid, plot_y_grid), method="nearest")
            grid[missing] = fill[missing]
    return plot_x, plot_y, plot_x_grid, plot_y_grid, u, v


def component_at_height_above_terrain(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    values: np.ndarray,
    terrain_interp: RegularGridInterpolator,
    lift: float,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
) -> np.ndarray:
    points, inverse = np.unique(
        np.column_stack((np.round(x, 8), np.round(y, 8))),
        axis=0,
        return_inverse=True,
    )
    interpolated = np.full(points.shape[0], np.nan)
    target_z = terrain_interp(points[:, [1, 0]]) + lift
    for index in range(points.shape[0]):
        keep = inverse == index
        z_column = z[keep]
        value_column = values[keep]
        order = np.argsort(z_column)
        z_unique, unique_indices = np.unique(z_column[order], return_index=True)
        value_unique = value_column[order][unique_indices]
        if target_z[index] < z_unique[0] or target_z[index] > z_unique[-1]:
            continue
        interpolated[index] = np.interp(target_z[index], z_unique, value_unique)

    valid = np.isfinite(interpolated)
    grid = griddata(points[valid], interpolated[valid], (x_grid, y_grid), method="linear")
    missing = np.isnan(grid)
    if np.any(missing):
        nearest = griddata(points[valid], interpolated[valid], (x_grid, y_grid), method="nearest")
        grid[missing] = nearest[missing]
    return grid


def trace_streamlines(
    seeds: np.ndarray,
    u_interp: RegularGridInterpolator,
    v_interp: RegularGridInterpolator,
    terrain_interp: RegularGridInterpolator,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    ground_lift: float,
    steps: int,
    step_length: float,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    paths: list[np.ndarray] = []
    path_speeds: list[np.ndarray] = []
    for seed in seeds:
        point = seed.copy()
        path = [point.copy()]
        speeds = [0.0]
        for _ in range(steps):
            where = np.array([[point[1], point[0]]])
            u = float(u_interp(where)[0])
            v = float(v_interp(where)[0])
            speed = float(np.hypot(u, v))
            if speed < 1.0e-8:
                break
            point[0] += step_length * u / speed
            point[1] += step_length * v / speed
            if point[0] < xlim[0] or point[0] > xlim[1] or point[1] < ylim[0] or point[1] > ylim[1]:
                break
            point[2] = float(terrain_interp([[point[1], point[0]]])[0] + ground_lift)
            path.append(point.copy())
            speeds.append(speed)
            if point[0] >= xlim[1] - 8.0:
                break
        if len(path) > 3:
            paths.append(np.asarray(path))
            path_speeds.append(np.asarray(speeds))
    return paths, path_speeds


def camera_basis(
    eye: np.ndarray,
    target: np.ndarray,
    up_hint: np.ndarray = np.array([0.0, 0.0, 1.0]),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    forward = target - eye
    forward /= np.linalg.norm(forward)
    right = np.cross(forward, up_hint)
    right /= np.linalg.norm(right)
    up = np.cross(right, forward)
    up /= np.linalg.norm(up)
    return right, up, forward


def rasterize_triangle(
    image: np.ndarray,
    depth: np.ndarray,
    triangle: np.ndarray,
    color: np.ndarray,
) -> None:
    xs = triangle[:, 0]
    ys = triangle[:, 1]
    min_x = max(int(np.floor(xs.min())), 0)
    max_x = min(int(np.ceil(xs.max())), image.shape[1] - 1)
    min_y = max(int(np.floor(ys.min())), 0)
    max_y = min(int(np.ceil(ys.max())), image.shape[0] - 1)
    if min_x > max_x or min_y > max_y:
        return

    x0, y0, z0 = triangle[0]
    x1, y1, z1 = triangle[1]
    x2, y2, z2 = triangle[2]
    denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
    if abs(denom) < 1.0e-10:
        return

    yy, xx = np.mgrid[min_y : max_y + 1, min_x : max_x + 1]
    px = xx + 0.5
    py = yy + 0.5
    w0 = ((y1 - y2) * (px - x2) + (x2 - x1) * (py - y2)) / denom
    w1 = ((y2 - y0) * (px - x2) + (x0 - x2) * (py - y2)) / denom
    w2 = 1.0 - w0 - w1
    mask = (w0 >= 0.0) & (w1 >= 0.0) & (w2 >= 0.0)
    if not np.any(mask):
        return

    tri_depth = w0 * z0 + w1 * z1 + w2 * z2
    target_depth = depth[min_y : max_y + 1, min_x : max_x + 1]
    visible = mask & (tri_depth < target_depth)
    if not np.any(visible):
        return
    target_image = image[min_y : max_y + 1, min_x : max_x + 1]
    target_image[visible] = color
    target_depth[visible] = tri_depth[visible]


def project_points(
    points: np.ndarray,
    eye: np.ndarray,
    right: np.ndarray,
    up: np.ndarray,
    forward: np.ndarray,
    x_min: float,
    y_min: float,
    scale: float,
    x_offset: float,
    y_offset: float,
    height: int,
    z_scale: float,
) -> np.ndarray:
    scaled_points = points.copy()
    scaled_points[:, 2] *= z_scale
    delta = scaled_points - eye
    camera = np.column_stack((delta @ right, delta @ up, delta @ forward))
    projected = np.empty_like(camera)
    projected[:, 0] = x_offset + scale * (camera[:, 0] - x_min)
    projected[:, 1] = height - (y_offset + scale * (camera[:, 1] - y_min))
    projected[:, 2] = camera[:, 2]
    return projected


def make_projection(
    faces: list[np.ndarray],
    width: int,
    height: int,
    elev: float,
    azim: float,
    z_scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float, float]:
    vertices = np.vstack([face_set.reshape(-1, 3) for face_set in faces])
    scaled_vertices = vertices.copy()
    scaled_vertices[:, 2] *= z_scale
    mins = vertices.min(axis=0)
    maxs = vertices.max(axis=0)
    scaled_mins = scaled_vertices.min(axis=0)
    scaled_maxs = scaled_vertices.max(axis=0)
    center = 0.5 * (scaled_mins + scaled_maxs)
    span = max(maxs[0] - mins[0], maxs[1] - mins[1])
    elev_rad = np.deg2rad(elev)
    azim_rad = np.deg2rad(azim)
    view_direction = np.array(
        [
            np.cos(elev_rad) * np.cos(azim_rad),
            np.cos(elev_rad) * np.sin(azim_rad),
            np.sin(elev_rad),
        ]
    )
    eye = center + span * view_direction
    target = center
    right, up, forward = camera_basis(eye, target)

    camera_xy = []
    for face_set in faces:
        vertices = face_set.reshape(-1, 3)
        vertices = vertices.copy()
        vertices[:, 2] *= z_scale
        delta = vertices - eye
        camera_xy.append(np.column_stack((delta @ right, delta @ up)))
    camera_xy = np.vstack(camera_xy)
    x_min, y_min = camera_xy.min(axis=0)
    x_max, y_max = camera_xy.max(axis=0)
    pad_x = 0.08 * (x_max - x_min)
    pad_y = 0.08 * (y_max - y_min)
    x_min -= pad_x
    x_max += pad_x
    y_min -= pad_y
    y_max += pad_y
    scale = min((width - 1) / (x_max - x_min), (height - 1) / (y_max - y_min))
    x_offset = 0.5 * (width - scale * (x_max - x_min))
    y_offset = 0.5 * (height - scale * (y_max - y_min))
    return eye, right, up, forward, x_min, y_min, scale, x_offset, y_offset


def shade_faces(faces: np.ndarray, base_colors: np.ndarray) -> np.ndarray:
    normals = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0])
    normal_lengths = np.linalg.norm(normals, axis=1)
    valid = normal_lengths > 1.0e-12
    normals[valid] /= normal_lengths[valid, None]
    light = np.array([-0.35, -0.5, 0.78])
    light /= np.linalg.norm(light)
    shade = 0.66 + 0.34 * np.abs(normals @ light)
    return np.clip(base_colors * shade[:, None], 0, 255).astype(np.uint8), valid


def draw_faces(
    image: np.ndarray,
    depth: np.ndarray,
    faces: np.ndarray,
    projected: np.ndarray,
    colors: np.ndarray,
    valid: np.ndarray,
) -> None:
    order = np.argsort(np.mean(projected[:, :, 2], axis=1))[::-1]
    for index in order:
        if valid[index]:
            rasterize_triangle(image, depth, projected[index], colors[index])


def main() -> None:
    terrain_faces = read_binary_stl(TERRAIN_FILE, max_faces=200000)
    building_faces = read_binary_stl(BUILDINGS_FILE, max_faces=200000)
    eye, right, up, forward, x_min, y_min, scale, x_offset, y_offset = make_projection(
        [terrain_faces, building_faces],
        WIDTH,
        HEIGHT,
        ELEVATION,
        AZIMUTH,
        Z_SCALE,
    )

    image = np.zeros((HEIGHT, WIDTH, 3), dtype=np.uint8)
    image[:, :] = np.array([238, 244, 247], dtype=np.uint8)
    depth = np.full((HEIGHT, WIDTH), np.inf, dtype=np.float64)

    terrain_z = np.mean(terrain_faces[:, :, 2], axis=1)
    terrain_norm = Normalize(vmin=float(np.nanmin(terrain_z)), vmax=float(np.nanmax(terrain_z)))
    terrain_base = LAND_CMAP(terrain_norm(terrain_z))[:, :3] * 255.0
    terrain_colors, terrain_valid = shade_faces(terrain_faces, terrain_base)
    terrain_projected = project_points(
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
        HEIGHT,
        Z_SCALE,
    ).reshape((-1, 3, 3))
    draw_faces(image, depth, terrain_faces, terrain_projected, terrain_colors, terrain_valid)

    building_base = np.full((len(building_faces), 3), 178.0)
    building_colors, building_valid = shade_faces(building_faces, building_base)
    building_projected = project_points(
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
        HEIGHT,
        Z_SCALE,
    ).reshape((-1, 3, 3))
    draw_faces(image, depth, building_faces, building_projected, building_colors, building_valid)

    x, y, z, u, v, _ = read_field(FIELD_FILE)
    terrain_x, terrain_y, terrain_grid_z = terrain_grid_from_stl(terrain_faces)
    terrain_interp = RegularGridInterpolator(
        (terrain_y[:, 0], terrain_x[0, :]),
        terrain_grid_z,
        bounds_error=False,
        fill_value=float(np.nanmin(terrain_grid_z)),
    )
    u_ground = component_at_height_above_terrain(
        x, y, z, u, terrain_interp, GROUND_LIFT, terrain_x, terrain_y
    )
    v_ground = component_at_height_above_terrain(
        x, y, z, v, terrain_interp, GROUND_LIFT, terrain_x, terrain_y
    )
    plot_x, plot_y, _, _, u_plot, v_plot = regularize_vector_field(
        terrain_x, terrain_y, u_ground, v_ground, n=180
    )
    u_interp = RegularGridInterpolator((plot_y, plot_x), u_plot, bounds_error=False, fill_value=0.0)
    v_interp = RegularGridInterpolator((plot_y, plot_x), v_plot, bounds_error=False, fill_value=0.0)
    xlim = (float(np.nanmin(terrain_x)), float(np.nanmax(terrain_x)))
    ylim = (float(np.nanmin(terrain_y)), float(np.nanmax(terrain_y)))
    source_x = xlim[0] + 12.0
    source_y = np.linspace(ylim[0] + 25.0, ylim[1] - 25.0, SOURCE_COUNT)
    seeds = np.column_stack(
        (
            np.full(SOURCE_COUNT, source_x),
            source_y,
            terrain_interp(np.column_stack((source_y, np.full(SOURCE_COUNT, source_x))))
            + GROUND_LIFT,
        )
    )
    paths, path_speeds = trace_streamlines(
        seeds,
        u_interp,
        v_interp,
        terrain_interp,
        xlim,
        ylim,
        GROUND_LIFT,
        STREAMLINE_STEPS,
        STREAMLINE_STEP,
    )

    draw = ImageDraw.Draw(Image.fromarray(image))
    speed_values = np.concatenate(path_speeds) if path_speeds else np.array([0.0, 1.0])
    speed_norm = Normalize(vmin=float(np.nanpercentile(speed_values, 5)), vmax=float(np.nanpercentile(speed_values, 95)))
    pil_image = Image.fromarray(image)
    draw = ImageDraw.Draw(pil_image)
    for path, speeds in zip(paths, path_speeds):
        projected = project_points(
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
            HEIGHT,
            Z_SCALE,
        )
        for index in range(len(projected) - 1):
            color = tuple((np.array(FLOW_HEAT_CMAP(speed_norm(speeds[min(index, len(speeds) - 1)]))[:3]) * 255).astype(np.uint8))
            draw.line(
                [
                    (float(projected[index, 0]), float(projected[index, 1])),
                    (float(projected[index + 1, 0]), float(projected[index + 1, 1])),
                ],
                fill=color,
                width=LINE_WIDTH,
            )

    pil_image.save(OUTPUT_FILE)
    print(f"Wrote {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
