#!/usr/bin/env python3
"""Build a coarse second-order hexahedral Gmsh mesh for Neko.

The mesh is terrain-following: the bottom boundary is sampled from
terrain_dem.xyz and the top is a flat plane. Horizontal cells are structured,
while the vertical spacing uses a stretching exponent to keep more elements near
the terrain/building boundary. The output is a Gmsh 2.2 ASCII mesh containing
HEX27 volume elements and QUAD9 boundary elements, which can be converted with
gmsh2nek and rea2nbin.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


ORIGIN_X = 676000.0
ORIGIN_Y = 6577000.0
PIPELINE_DIR = Path(__file__).resolve().parents[1]
DEM_FILE = PIPELINE_DIR / "generated" / "terrain_dem.xyz"
OUTPUT_FILE = PIPELINE_DIR / "generated" / "urban_terrain.msh"
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


def make_coordinates(
    dem_x: np.ndarray,
    dem_y: np.ndarray,
    dem_z: np.ndarray,
    nx: int,
    ny: int,
    nz: int,
    top_z: float,
    stretch: float,
) -> tuple[np.ndarray, dict[tuple[int, int, int], int]]:
    x_coords = np.linspace(float(dem_x.min()), float(dem_x.max()), 2 * nx + 1)
    y_coords = np.linspace(float(dem_y.min()), float(dem_y.max()), 2 * ny + 1)
    eta = np.linspace(0.0, 1.0, 2 * nz + 1) ** stretch

    nodes = []
    node_ids: dict[tuple[int, int, int], int] = {}
    for k, s in enumerate(eta):
        for j, y in enumerate(y_coords):
            for i, x in enumerate(x_coords):
                terrain_z = interp_terrain(x, y, dem_x, dem_y, dem_z)
                z = terrain_z + s * (top_z - terrain_z)
                node_ids[(i, j, k)] = len(nodes) + 1
                nodes.append((x, y, z))
    return np.asarray(nodes), node_ids


def nid(node_ids: dict[tuple[int, int, int], int], i: int, j: int, k: int) -> int:
    return node_ids[(i, j, k)]


def hex27_nodes(
    node_ids: dict[tuple[int, int, int], int], ex: int, ey: int, ez: int
) -> list[int]:
    i = 2 * ex
    j = 2 * ey
    k = 2 * ez
    return [
        nid(node_ids, i, j, k),
        nid(node_ids, i + 2, j, k),
        nid(node_ids, i + 2, j + 2, k),
        nid(node_ids, i, j + 2, k),
        nid(node_ids, i, j, k + 2),
        nid(node_ids, i + 2, j, k + 2),
        nid(node_ids, i + 2, j + 2, k + 2),
        nid(node_ids, i, j + 2, k + 2),
        nid(node_ids, i + 1, j, k),
        nid(node_ids, i, j + 1, k),
        nid(node_ids, i, j, k + 1),
        nid(node_ids, i + 2, j + 1, k),
        nid(node_ids, i + 2, j, k + 1),
        nid(node_ids, i + 1, j + 2, k),
        nid(node_ids, i + 2, j + 2, k + 1),
        nid(node_ids, i, j + 2, k + 1),
        nid(node_ids, i + 1, j, k + 2),
        nid(node_ids, i, j + 1, k + 2),
        nid(node_ids, i + 2, j + 1, k + 2),
        nid(node_ids, i + 1, j + 2, k + 2),
        nid(node_ids, i + 1, j + 1, k),
        nid(node_ids, i + 1, j, k + 1),
        nid(node_ids, i, j + 1, k + 1),
        nid(node_ids, i + 2, j + 1, k + 1),
        nid(node_ids, i + 1, j + 2, k + 1),
        nid(node_ids, i + 1, j + 1, k + 2),
        nid(node_ids, i + 1, j + 1, k + 1),
    ]


def quad9(
    node_ids: dict[tuple[int, int, int], int],
    corners: tuple[tuple[int, int, int], ...],
    edges: tuple[tuple[int, int, int], ...],
    center: tuple[int, int, int],
) -> list[int]:
    return [nid(node_ids, *p) for p in (*corners, *edges, center)]


def boundary_elements(
    node_ids: dict[tuple[int, int, int], int], nx: int, ny: int, nz: int
) -> list[tuple[int, list[int]]]:
    elements: list[tuple[int, list[int]]] = []

    for ey in range(ny):
        for ez in range(nz):
            j = 2 * ey
            k = 2 * ez
            elements.append(
                (
                    1,
                    quad9(
                        node_ids,
                        ((0, j, k), (0, j, k + 2), (0, j + 2, k + 2), (0, j + 2, k)),
                        ((0, j, k + 1), (0, j + 1, k + 2), (0, j + 2, k + 1), (0, j + 1, k)),
                        (0, j + 1, k + 1),
                    ),
                )
            )
            i = 2 * nx
            elements.append(
                (
                    2,
                    quad9(
                        node_ids,
                        ((i, j, k), (i, j + 2, k), (i, j + 2, k + 2), (i, j, k + 2)),
                        ((i, j + 1, k), (i, j + 2, k + 1), (i, j + 1, k + 2), (i, j, k + 1)),
                        (i, j + 1, k + 1),
                    ),
                )
            )

    for ex in range(nx):
        for ez in range(nz):
            i = 2 * ex
            k = 2 * ez
            elements.append(
                (
                    3,
                    quad9(
                        node_ids,
                        ((i, 0, k), (i + 2, 0, k), (i + 2, 0, k + 2), (i, 0, k + 2)),
                        ((i + 1, 0, k), (i + 2, 0, k + 1), (i + 1, 0, k + 2), (i, 0, k + 1)),
                        (i + 1, 0, k + 1),
                    ),
                )
            )
            j = 2 * ny
            elements.append(
                (
                    4,
                    quad9(
                        node_ids,
                        ((i, j, k), (i, j, k + 2), (i + 2, j, k + 2), (i + 2, j, k)),
                        ((i, j, k + 1), (i + 1, j, k + 2), (i + 2, j, k + 1), (i + 1, j, k)),
                        (i + 1, j, k + 1),
                    ),
                )
            )

    for ex in range(nx):
        for ey in range(ny):
            i = 2 * ex
            j = 2 * ey
            elements.append(
                (
                    5,
                    quad9(
                        node_ids,
                        ((i, j, 0), (i, j + 2, 0), (i + 2, j + 2, 0), (i + 2, j, 0)),
                        ((i, j + 1, 0), (i + 1, j + 2, 0), (i + 2, j + 1, 0), (i + 1, j, 0)),
                        (i + 1, j + 1, 0),
                    ),
                )
            )
            k = 2 * nz
            elements.append(
                (
                    6,
                    quad9(
                        node_ids,
                        ((i, j, k), (i + 2, j, k), (i + 2, j + 2, k), (i, j + 2, k)),
                        ((i + 1, j, k), (i + 2, j + 1, k), (i + 1, j + 2, k), (i, j + 1, k)),
                        (i + 1, j + 1, k),
                    ),
                )
            )
    return elements


def write_msh(
    path: Path,
    nodes: np.ndarray,
    node_ids: dict[tuple[int, int, int], int],
    nx: int,
    ny: int,
    nz: int,
) -> None:
    boundary = boundary_elements(node_ids, nx, ny, nz)
    n_hex = nx * ny * nz
    total_elements = len(boundary) + n_hex

    with path.open("w", encoding="ascii") as handle:
        handle.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        handle.write("$PhysicalNames\n7\n")
        for tag, name in PHYSICAL_NAMES.items():
            handle.write(f'2 {tag} "{name}"\n')
        handle.write('3 1 "fluid"\n')
        handle.write("$EndPhysicalNames\n")

        handle.write(f"$Nodes\n{len(nodes)}\n")
        for idx, (x, y, z) in enumerate(nodes, start=1):
            handle.write(f"{idx} {x:.12g} {y:.12g} {z:.12g}\n")
        handle.write("$EndNodes\n")

        handle.write(f"$Elements\n{total_elements}\n")
        element_id = 1
        for physical_id, elem_nodes in boundary:
            handle.write(
                f"{element_id} 10 2 {physical_id} {physical_id} "
                + " ".join(str(n) for n in elem_nodes)
                + "\n"
            )
            element_id += 1

        for ez in range(nz):
            for ey in range(ny):
                for ex in range(nx):
                    elem_nodes = hex27_nodes(node_ids, ex, ey, ez)
                    handle.write(
                        f"{element_id} 12 2 1 1 "
                        + " ".join(str(n) for n in elem_nodes)
                        + "\n"
                    )
                    element_id += 1
        handle.write("$EndElements\n")


def main() -> None:
    dem_x, dem_y, dem_z = read_dem(DEM_FILE)

    nodes, node_ids = make_coordinates(
        dem_x,
        dem_y,
        dem_z,
        NX,
        NY,
        NZ,
        TOP_Z,
        VERTICAL_STRETCH,
    )
    write_msh(OUTPUT_FILE, nodes, node_ids, NX, NY, NZ)
    print(
        f"Wrote {OUTPUT_FILE}: {NX * NY * NZ} HEX27 elements, "
        f"{len(nodes)} nodes"
    )


if __name__ == "__main__":
    main()
