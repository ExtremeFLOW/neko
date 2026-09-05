#!/usr/bin/env python3
"""Create a coordinate-only Neko field from an nmsh mesh.

This is useful for Brinkman-cache generation when we need the physical GLL
coordinates but do not want to run a solver just to write an initial field.
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np


HEADER_SIZE = 132
MARKER = 6.54321


def read_nmsh(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open("rb") as handle:
        n_elements, gdim = struct.unpack("<2i", handle.read(8))
        if gdim != 3:
            raise RuntimeError(f"Expected 3D nmsh, got gdim={gdim}")
        element_ids = np.empty(n_elements, dtype=np.int32)
        vertices = np.empty((n_elements, 8, 3), dtype=np.float64)
        for e in range(n_elements):
            element_ids[e] = struct.unpack("<i", handle.read(4))[0]
            for v in range(8):
                _, x, y, z = struct.unpack("<i3d", handle.read(28))
                vertices[e, v, :] = (x, y, z)
    return element_ids, vertices


def gll_nodes(order: int) -> np.ndarray:
    if order < 1:
        raise RuntimeError("Polynomial order must be at least 1.")
    if order == 1:
        return np.array([-1.0, 1.0], dtype=np.float64)
    poly = np.polynomial.legendre.Legendre.basis(order)
    interior = poly.deriv().roots()
    return np.concatenate(([-1.0], np.sort(interior), [1.0])).astype(np.float64)


def trilinear_hex_points(vertices: np.ndarray, order: int) -> np.ndarray:
    nodes = gll_nodes(order)
    n = nodes.size
    points = np.empty((vertices.shape[0], n * n * n, 3), dtype=np.float64)
    signs = np.array(
        [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    idx = 0
    for zeta in nodes:
        for eta in nodes:
            for xi in nodes:
                weights = (
                    0.125
                    * (1.0 + signs[:, 0] * xi)
                    * (1.0 + signs[:, 1] * eta)
                    * (1.0 + signs[:, 2] * zeta)
                )
                points[:, idx, :] = weights @ vertices
                idx += 1
    return points


def write_field(path: Path, element_ids: np.ndarray, points: np.ndarray, order: int) -> None:
    lx = order + 1
    ly = lx
    lz = lx
    nelv = int(points.shape[0])
    nelgt = nelv
    npts = lx * ly * lz
    header_text = (
        f"#std 8 {lx:2d} {ly:2d} {lz:2d} {nelv:10d} {nelgt:10d} "
        f"{0.0:20.13E} {0:9d} {1:6d} {1:6d} {'X':<10s}"
    )
    header_bytes = header_text.encode("ascii")
    if len(header_bytes) > HEADER_SIZE:
        raise RuntimeError("Generated Neko field header is too long")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        handle.write(header_bytes.ljust(HEADER_SIZE, b" "))
        handle.write(struct.pack("<f", MARKER))
        element_ids.astype("<i4", copy=False).tofile(handle)
        flat = points.reshape(nelv, npts, 3)
        for element in range(nelv):
            flat[element, :, 0].astype("<f8", copy=False).tofile(handle)
            flat[element, :, 1].astype("<f8", copy=False).tofile(handle)
            flat[element, :, 2].astype("<f8", copy=False).tofile(handle)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--polynomial-order", type=int, default=2)
    args = parser.parse_args()

    element_ids, vertices = read_nmsh(args.mesh)
    points = trilinear_hex_points(vertices, args.polynomial_order)
    write_field(args.output, element_ids, points, args.polynomial_order)
    print(f"Wrote {args.output} with {len(element_ids)} elements and {(args.polynomial_order + 1) ** 3} points/element")


if __name__ == "__main__":
    main()
