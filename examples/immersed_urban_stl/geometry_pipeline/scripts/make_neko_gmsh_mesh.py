#!/usr/bin/env python3
"""Build a terrain-following hexahedral Neko mesh.

The mesh is terrain-following: the bottom boundary is sampled from
terrain_dem.xyz and the top is a flat plane. Horizontal cells are structured,
while the vertical spacing uses a stretching exponent to keep more elements near
the terrain/building boundary. The output is a native Neko .nmsh file.
"""

from __future__ import annotations

from pathlib import Path
import struct

import numpy as np


ORIGIN_X = 676000.0
ORIGIN_Y = 6577000.0
PIPELINE_DIR = Path(__file__).resolve().parents[1]
DEM_FILE = PIPELINE_DIR / "generated" / "terrain_dem.xyz"
OUTPUT_FILE = PIPELINE_DIR / "generated" / "urban_terrain.nmsh"
NX = 144
NY = 144
NZ = 45
TOP_Z = 250.0
VERTICAL_STRETCH = 1.8

PHYSICAL_NAMES = {
    1: "inlet_xmin",
    2: "outlet_xmax",
    3: "side_ymin",
    4: "side_ymax",
    5: "terrain",
    6: "top",
}


def read_dem(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    x = data[:, 0] - ORIGIN_X
    y = data[:, 1] - ORIGIN_Y
    z = data[:, 2]

    xs = np.unique(np.round(x, 10))
    ys = np.unique(np.round(y, 10))
    nx = len(xs)
    ny = len(ys)

    order = np.lexsort((x, y))
    sorted_x = np.round(x[order], 10)
    sorted_y = np.round(y[order], 10)
    sorted_z = z[order]

    grid = np.empty((ny, nx))
    for idx, (_, yy) in enumerate(zip(sorted_x, sorted_y)):
        row = np.where(ys == yy)[0][0]
        col = idx % nx
        grid[row, col] = sorted_z[idx]

    return xs, ys, grid


def interp_terrain(
    x: float, y: float, dem_x: np.ndarray, dem_y: np.ndarray, dem_z: np.ndarray
) -> float:
    i = int(np.searchsorted(dem_x, x) - 1)
    j = int(np.searchsorted(dem_y, y) - 1)
    i = max(0, min(i, len(dem_x) - 2))
    j = max(0, min(j, len(dem_y) - 2))

    x0 = dem_x[i]
    x1 = dem_x[i + 1]
    y0 = dem_y[j]
    y1 = dem_y[j + 1]
    tx = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
    ty = 0.0 if y1 == y0 else (y - y0) / (y1 - y0)

    z00 = dem_z[j, i]
    z10 = dem_z[j, i + 1]
    z01 = dem_z[j + 1, i]
    z11 = dem_z[j + 1, i + 1]
    return float(
        (1.0 - tx) * (1.0 - ty) * z00
        + tx * (1.0 - ty) * z10
        + (1.0 - tx) * ty * z01
        + tx * ty * z11
    )


def vertex_id(i: int, j: int, k: int, nx: int, ny: int) -> int:
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


def element_id(ex: int, ey: int, ez: int, nx: int, ny: int) -> int:
    return 1 + ex + nx * (ey + ny * ez)


def write_zone(handle, element: int, face: int, label: int) -> None:
    handle.write(struct.pack("<9i", element, face, 0, label, 0, 0, 0, 0, 7))


def write_nmsh(
    path: Path,
    dem_x: np.ndarray,
    dem_y: np.ndarray,
    dem_z: np.ndarray,
    nx: int,
    ny: int,
    nz: int,
    top_z: float,
    stretch: float,
) -> None:
    x_coords = np.linspace(float(dem_x.min()), float(dem_x.max()), nx + 1)
    y_coords = np.linspace(float(dem_y.min()), float(dem_y.max()), ny + 1)
    eta = np.linspace(0.0, 1.0, nz + 1) ** stretch

    terrain = np.empty((ny + 1, nx + 1))
    for j, y in enumerate(y_coords):
        for i, x in enumerate(x_coords):
            terrain[j, i] = interp_terrain(x, y, dem_x, dem_y, dem_z)

    def vertex(i: int, j: int, k: int) -> tuple[int, float, float, float]:
        z0 = terrain[j, i]
        z = z0 + eta[k] * (top_z - z0)
        return vertex_id(i, j, k, nx, ny), x_coords[i], y_coords[j], z

    n_elements = nx * ny * nz
    n_zones = 2 * (ny * nz + nx * nz + nx * ny)

    with path.open("wb") as handle:
        handle.write(struct.pack("<2i", n_elements, 3))

        for ez in range(nz):
            for ey in range(ny):
                for ex in range(nx):
                    handle.write(struct.pack("<i", element_id(ex, ey, ez, nx, ny)))
                    for idx, x, y, z in (
                        vertex(ex, ey, ez),
                        vertex(ex + 1, ey, ez),
                        vertex(ex + 1, ey + 1, ez),
                        vertex(ex, ey + 1, ez),
                        vertex(ex, ey, ez + 1),
                        vertex(ex + 1, ey, ez + 1),
                        vertex(ex + 1, ey + 1, ez + 1),
                        vertex(ex, ey + 1, ez + 1),
                    ):
                        handle.write(struct.pack("<i3d", idx, x, y, z))

        handle.write(struct.pack("<i", n_zones))

        for ey in range(ny):
            for ez in range(nz):
                write_zone(handle, element_id(0, ey, ez, nx, ny), 1, 1)
                write_zone(handle, element_id(nx - 1, ey, ez, nx, ny), 2, 2)

        for ex in range(nx):
            for ez in range(nz):
                write_zone(handle, element_id(ex, 0, ez, nx, ny), 3, 3)
                write_zone(handle, element_id(ex, ny - 1, ez, nx, ny), 4, 4)

        for ex in range(nx):
            for ey in range(ny):
                write_zone(handle, element_id(ex, ey, 0, nx, ny), 5, 5)
                write_zone(handle, element_id(ex, ey, nz - 1, nx, ny), 6, 6)

        handle.write(struct.pack("<i", 0))


def main() -> None:
    dem_x, dem_y, dem_z = read_dem(DEM_FILE)

    write_nmsh(OUTPUT_FILE, dem_x, dem_y, dem_z, NX, NY, NZ, TOP_Z, VERTICAL_STRETCH)
    print(
        f"Wrote {OUTPUT_FILE}: {NX * NY * NZ} hex elements, "
        f"{(NX + 1) * (NY + 1) * (NZ + 1)} vertices"
    )


if __name__ == "__main__":
    main()
