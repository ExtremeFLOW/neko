#!/usr/bin/env python3
"""Check actual cylindrical side labels in a little-endian Neko mesh."""

import argparse
import json
import math
from pathlib import Path
import struct

import numpy as np


def check_boundary_zones(path: Path, inflow_from: float, inflow_width: float,
                         outflow_width: float) -> dict:
    path = Path(path)
    if not 0.0 < inflow_width < 360.0 or not 0.0 < outflow_width < 360.0:
        raise ValueError("Boundary arc widths must be between 0 and 360 degrees")
    vertex = np.dtype([("id", "<i4"), ("xyz", "<f8", (3,))])
    element = np.dtype([("id", "<i4"), ("vertices", vertex, (8,))])
    with path.open("rb") as f:
        count, dim = struct.unpack("<2i", f.read(8))
        if dim != 3 or count <= 0:
            raise ValueError("Expected a nonempty 3D little-endian .nmsh")
        f.seek(8 + count * element.itemsize)
        nzone, = struct.unpack("<i", f.read(4))
        zones = np.fromfile(f, dtype="<i4", count=nzone * 9).reshape(nzone, 9)
    elements = np.memmap(path, mode="r", offset=8, dtype=element, shape=(count,))
    side_edges = {1: (0, 3), 2: (1, 2), 3: (0, 1), 4: (3, 2)}
    edges = {}
    for row in zones:
        eid, face, label = int(row[0]), int(row[1]), int(row[3])
        if face not in side_edges:
            continue
        if not 1 <= eid <= count or label not in (1, 2, 3):
            raise ValueError(f"Invalid cylindrical side label: element={eid}, zone={label}")
        if elements[eid - 1]["id"] != eid:
            raise ValueError("Expected sequential element IDs from the cylinder generator")
        xy = elements[eid - 1]["vertices"]["xyz"][list(side_edges[face]), :2]
        key = tuple(sorted(tuple(np.round(p, 8)) for p in xy))
        if key in edges and edges[key] != label:
            raise ValueError("Side labels change between vertical layers")
        edges[key] = label

    widths = {1: 0.0, 2: 0.0, 3: 0.0}
    inlet_vertices = {}
    transitions = [(inflow_from + sign * inflow_width / 2.0) % 360.0 for sign in (-1, 1)]
    transitions += [(inflow_from + 180.0 + sign * outflow_width / 2.0) % 360.0 for sign in (-1, 1)]
    distance = lambda a, b: abs((a - b + 180.0) % 360.0 - 180.0)
    for (a, b), label in edges.items():
        ba, bb = (math.degrees(math.atan2(p[0], p[1])) % 360.0 for p in (a, b))
        delta = (bb - ba + 180.0) % 360.0 - 180.0
        mid = (ba + delta / 2.0) % 360.0
        expected = (1 if distance(mid, inflow_from) <= inflow_width / 2.0 else
                    2 if distance(mid, inflow_from + 180.0) <= outflow_width / 2.0 else 3)
        if label != expected:
            raise ValueError(f"Side zone mismatch at bearing {mid:.6f}: got {label}, expected {expected}")
        if any(distance(mid, t) < abs(delta) / 2.0 - 1e-6 for t in transitions):
            raise ValueError(f"Boundary face crosses a requested zone endpoint near {mid:.6f} degrees")
        widths[label] += abs(delta)
        if label == 1:
            for p in (a, b):
                inlet_vertices[p] = inlet_vertices.get(p, 0) + 1
    if not math.isclose(sum(widths.values()), 360.0, abs_tol=1e-5):
        raise ValueError("Cylindrical boundary does not span 360 degrees")
    if not math.isclose(widths[1], inflow_width, abs_tol=1e-5):
        raise ValueError(f"Actual inlet {widths[1]:.6f} degrees differs from requested {inflow_width}")
    endpoints = sorted(math.degrees(math.atan2(p[0], p[1])) % 360.0
                       for p, degree in inlet_vertices.items() if degree == 1)
    expected_endpoints = transitions[:2]
    if len(endpoints) != 2 or any(min(distance(t, p) for p in endpoints) > 1e-5 for t in expected_endpoints):
        raise ValueError(f"Incorrect inlet endpoints: {endpoints}")
    return {"inlet_width_degrees": widths[1], "inlet_endpoints_degrees": endpoints,
            "side_zone_widths_degrees": widths, "horizontal_boundary_edges": len(edges)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--inflow-from", type=float, default=225.0)
    parser.add_argument("--inflow-width", type=float, required=True)
    parser.add_argument("--outflow-width", type=float, default=90.0)
    args = parser.parse_args()
    print(json.dumps(check_boundary_zones(args.mesh, args.inflow_from, args.inflow_width,
                                         args.outflow_width), indent=2))
