#!/usr/bin/env python3
"""Build a local Södermalm cylindrical hex-mesh prototype.

This version uses Gmsh to create a proper quad mesh of a circular disk, then
extrudes those quads into Neko hex elements. The terrain is real/interpolated
over Södermalm and smoothly relaxes to the inferred water level outside the
island, so the circular boundary sits at water level.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import struct
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT_GPKG = HERE / "input" / "inner_stockholm_15km_layers_clipped.gpkg"
OUT = HERE / "generated"
FIG = HERE / "figures"

OGR2OGR = os.environ.get("OGR2OGR", shutil.which("ogr2ogr") or "/opt/homebrew/bin/ogr2ogr")
GMSH = os.environ.get("GMSH", shutil.which("gmsh") or "/opt/homebrew/bin/gmsh")
RSVG = os.environ.get("RSVG", shutil.which("rsvg-convert") or "/opt/homebrew/bin/rsvg-convert")
OSGEO_PYTHON = os.environ.get("OSGEO_PYTHON", shutil.which("python3") or "/opt/homebrew/bin/python3")
REPO_MESH_CHECKER = HERE.parents[2] / ".codex_build" / "cpu_install" / "bin" / "mesh_checker"
MESH_CHECKER = str(REPO_MESH_CHECKER) if REPO_MESH_CHECKER.exists() else (shutil.which("mesh_checker") or "/usr/local/bin/mesh_checker")

NZ = 10
BUFFER_M = 220.0
DOMAIN_HEIGHT_ABOVE_LOWEST_M = 1000.0
VERTICAL_STRETCH = 1.5
MESH_SIZE_M = 260.0
LAND_MESH_SIZE_M = 100.0
WATER_MESH_SIZE_M = 280.0
SHORELINE_SIMPLIFY_M = 20.0
SHORE_BLEND_M = 180.0
INFLOW_FROM_DEG = 225.0
INFLOW_ARC_WIDTH_DEG = 180.0
OUTFLOW_ARC_WIDTH_DEG = 90.0
CIRCLE_ARC_STEP_DEG = 45.0
BUILDING_MIN_TERRAIN_ABOVE_WATER_M = 0.75
BUILDING_BASE_CLEARANCE_M = 0.5
BUILDING_MIN_FOOTPRINT_AREA_M2 = 0.0
BUILDING_MIN_EDGE_M = 0.0
BUILDING_SIMPLIFY_M = 0.0
BUILDING_MIN_HEIGHT_M = 2.0
BUILDING_MAX_HEIGHT_M = 65.0


# Rough hint used by the mask-derivation helper to find the correct connected
# land component. It is not used as the final shoreline.
SODERMALM_WGS84 = [
    (17.9985, 59.3150),
    (18.0120, 59.3215),
    (18.0400, 59.3228),
    (18.0730, 59.3215),
    (18.1015, 59.3130),
    (18.1065, 59.3035),
    (18.0890, 59.2965),
    (18.0630, 59.2950),
    (18.0330, 59.3010),
    (18.0100, 59.3065),
    (17.9985, 59.3150),
]


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_mesh_checker(mesh_path: Path) -> str:
    env = os.environ.copy()
    lib_paths = ["/usr/local/neko_install/bin", "/usr/local/json-fortran_install/lib"]
    env["DYLD_LIBRARY_PATH"] = ":".join(lib_paths + ([env["DYLD_LIBRARY_PATH"]] if env.get("DYLD_LIBRARY_PATH") else []))
    print("+", f"(cd {mesh_path.parent} && {MESH_CHECKER} {mesh_path.name})")
    result = subprocess.run(
        [MESH_CHECKER, mesh_path.name],
        check=True,
        env=env,
        cwd=mesh_path.parent,
        text=True,
        capture_output=True,
    )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    return result.stdout + result.stderr


def load_geojson(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_geojson(path: Path, geometry: dict, crs_name: str = "EPSG:4326") -> None:
    data = {
        "type": "FeatureCollection",
        "name": path.stem,
        "crs": {"type": "name", "properties": {"name": crs_name}},
        "features": [{"type": "Feature", "properties": {}, "geometry": geometry}],
    }
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def polygon_rings(geometry: dict) -> list[list[tuple[float, float]]]:
    if not geometry:
        return []
    if geometry["type"] == "Polygon":
        return [[(float(x), float(y)) for x, y, *_ in geometry["coordinates"][0]]]
    if geometry["type"] == "MultiPolygon":
        return [[(float(x), float(y)) for x, y, *_ in poly[0]] for poly in geometry["coordinates"]]
    return []


def polygon_boundary_rings(geometry: dict) -> list[list[tuple[float, float]]]:
    if not geometry:
        return []
    if geometry["type"] == "Polygon":
        return [[(float(x), float(y)) for x, y, *_ in ring] for ring in geometry["coordinates"]]
    if geometry["type"] == "MultiPolygon":
        return [[(float(x), float(y)) for x, y, *_ in ring] for poly in geometry["coordinates"] for ring in poly]
    return []


def line_points(geometry: dict) -> list[tuple[float, float]]:
    if not geometry:
        return []
    if geometry["type"] == "LineString":
        return [(float(x), float(y)) for x, y, *_ in geometry["coordinates"]]
    if geometry["type"] == "MultiLineString":
        return [(float(x), float(y)) for part in geometry["coordinates"] for x, y, *_ in part]
    return []


def point_in_polygon(x: float, y: float, ring: list[tuple[float, float]]) -> bool:
    inside = False
    for i, (x0, y0) in enumerate(ring):
        x1, y1 = ring[(i + 1) % len(ring)]
        if (y0 > y) != (y1 > y):
            x_cross = (x1 - x0) * (y - y0) / (y1 - y0 + 1e-300) + x0
            if x < x_cross:
                inside = not inside
    return inside


def point_in_land_mask(x: float, y: float, boundary_rings: list[list[tuple[float, float]]]) -> bool:
    # Odd-even containment handles both polygon holes and disjoint island rings.
    inside = False
    for ring in boundary_rings:
        if point_in_polygon(x, y, ring):
            inside = not inside
    return inside


def point_segment_distance(px: float, py: float, ax: float, ay: float, bx: float, by: float) -> float:
    vx, vy = bx - ax, by - ay
    wx, wy = px - ax, py - ay
    denom = vx * vx + vy * vy
    t = 0.0 if denom == 0.0 else max(0.0, min(1.0, (wx * vx + wy * vy) / denom))
    qx, qy = ax + t * vx, ay + t * vy
    return math.hypot(px - qx, py - qy)


def distance_to_ring(x: float, y: float, ring: list[tuple[float, float]]) -> float:
    return min(point_segment_distance(x, y, *ring[i], *ring[(i + 1) % len(ring)]) for i in range(len(ring) - 1))


def smoothstep(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


def terrain_samples(contours: Path, center: tuple[float, float]) -> tuple[np.ndarray, float]:
    rows: list[tuple[float, float, float]] = []
    for feature in load_geojson(contours)["features"]:
        z = feature["properties"].get("HOJD")
        if z is None:
            continue
        for x, y in line_points(feature["geometry"]):
            rows.append((x - center[0], y - center[1], float(z)))
    if len(rows) < 10:
        raise SystemExit("Too few terrain contour samples for Södermalm cylinder.")
    data = np.array(rows, dtype=float)
    water_level = float(np.min(data[:, 2]))
    return data, water_level


def idw_height(x: float, y: float, samples: np.ndarray) -> float:
    xy = samples[:, :2]
    zz = samples[:, 2]
    dist2 = (xy[:, 0] - x) ** 2 + (xy[:, 1] - y) ** 2
    idx = np.argpartition(dist2, min(20, len(dist2) - 1))[:20]
    d = np.sqrt(dist2[idx])
    if float(d.min()) < 1e-9:
        return float(zz[idx[int(d.argmin())]])
    w = 1.0 / np.maximum(d, 2.0) ** 2
    return float(np.sum(w * zz[idx]) / np.sum(w))


def terrain_height(
    x: float,
    y: float,
    samples: np.ndarray,
    water_level: float,
    shoreline: list[list[tuple[float, float]]],
    center: tuple[float, float],
    radius: float,
) -> float:
    gx, gy = x + center[0], y + center[1]
    inside = point_in_land_mask(gx, gy, shoreline)
    if not inside:
        return water_level
    d = min(distance_to_ring(gx, gy, ring) for ring in shoreline)
    weight = smoothstep(d / SHORE_BLEND_M)
    raw = max(water_level, idw_height(x, y, samples))
    # Keep the outer cylinder exactly flat water even if the mask was too large.
    edge = math.hypot(x, y)
    edge_weight = smoothstep(max(0.0, radius - edge) / max(SHORE_BLEND_M, 1.0))
    weight *= edge_weight
    return water_level + weight * (raw - water_level)


def perpendicular_distance(point: tuple[float, float], start: tuple[float, float], end: tuple[float, float]) -> float:
    return point_segment_distance(point[0], point[1], start[0], start[1], end[0], end[1])


def simplify_ring(points: list[tuple[float, float]], tolerance: float) -> list[tuple[float, float]]:
    pts = points[:-1] if points and points[0] == points[-1] else list(points)
    if len(pts) <= 4:
        return pts

    def rdp(rows: list[tuple[float, float]]) -> list[tuple[float, float]]:
        if len(rows) <= 2:
            return rows
        start, end = rows[0], rows[-1]
        distances = [perpendicular_distance(row, start, end) for row in rows[1:-1]]
        if not distances:
            return rows
        idx = int(np.argmax(distances)) + 1
        if distances[idx - 1] <= tolerance:
            return [start, end]
        return rdp(rows[: idx + 1])[:-1] + rdp(rows[idx:])

    simplified = rdp(pts + [pts[0]])[:-1]
    if len(simplified) < 4:
        return pts
    return simplified


def compass_bearing_deg(x: float, y: float) -> float:
    return (math.degrees(math.atan2(x, y)) + 360.0) % 360.0


def circular_distance_deg(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def is_inflow_bearing(bearing: float) -> bool:
    return circular_distance_deg(bearing, INFLOW_FROM_DEG) <= 0.5 * INFLOW_ARC_WIDTH_DEG + 1e-9


def is_outflow_bearing(bearing: float) -> bool:
    flow_to_deg = (INFLOW_FROM_DEG + 180.0) % 360.0
    return circular_distance_deg(bearing, flow_to_deg) <= 0.5 * OUTFLOW_ARC_WIDTH_DEG + 1e-9


def cylindrical_side_zone(bearing: float) -> int:
    if is_inflow_bearing(bearing):
        return 1
    if is_outflow_bearing(bearing):
        return 2
    return 3


def circle_point(angle_deg: float, radius: float) -> tuple[float, float]:
    angle = math.radians(angle_deg)
    return radius * math.cos(angle), radius * math.sin(angle)


def inlet_arc_polylines(radius: float, n: int = 361) -> list[list[tuple[float, float]]]:
    angles = np.linspace(0.0, 360.0, n, endpoint=False)
    runs: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []
    for angle in angles:
        x, y = circle_point(float(angle), radius)
        if is_inflow_bearing(compass_bearing_deg(x, y)):
            current.append((x, y))
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    if len(runs) > 1 and is_inflow_bearing(compass_bearing_deg(*circle_point(0.0, radius))):
        runs[0] = runs[-1] + runs[0]
        runs.pop()
    return runs


def write_gmsh_geo(path: Path, radius: float, shoreline: list[list[tuple[float, float]]], center: tuple[float, float]) -> None:
    if 360.0 % CIRCLE_ARC_STEP_DEG != 0.0:
        raise ValueError("CIRCLE_ARC_STEP_DEG must divide 360 degrees exactly")
    n_arcs = int(round(360.0 / CIRCLE_ARC_STEP_DEG))
    circle_points = [circle_point(i * CIRCLE_ARC_STEP_DEG, radius) for i in range(n_arcs)]
    side_curves: dict[int, list[int]] = {1: [], 2: [], 3: []}

    lines: list[str] = [
        'SetFactory("OpenCASCADE");',
        "Mesh.MshFileVersion = 2.2;",
        "Mesh.Algorithm = 8;",
        "Mesh.RecombineAll = 1;",
        "Mesh.SubdivisionAlgorithm = 1;",
        "Mesh.Optimize = 1;",
        f"Mesh.CharacteristicLengthMin = {LAND_MESH_SIZE_M};",
        f"Mesh.CharacteristicLengthMax = {WATER_MESH_SIZE_M};",
        f"Point(1) = {{0, 0, 0, {WATER_MESH_SIZE_M}}};",
    ]
    for i, (x, y) in enumerate(circle_points, start=2):
        lines.append(f"Point({i}) = {{{x}, {y}, 0, {WATER_MESH_SIZE_M}}};")
    for i in range(n_arcs):
        curve_id = i + 1
        start = i + 2
        end = 2 if i == n_arcs - 1 else i + 3
        lines.append(f"Circle({curve_id}) = {{{start}, 1, {end}}};")
        mid_angle = (i + 0.5) * CIRCLE_ARC_STEP_DEG
        mx, my = circle_point(mid_angle, 1.0)
        side_curves[cylindrical_side_zone(compass_bearing_deg(mx, my))].append(curve_id)
    lines.append(f"Curve Loop(1) = {{{', '.join(str(i + 1) for i in range(n_arcs))}}};")

    next_point = 100
    next_curve = 100
    next_loop = 10

    def add_ring(ring: list[tuple[float, float]], lc: float) -> int:
        nonlocal next_point, next_curve, next_loop
        simple = simplify_ring(ring, SHORELINE_SIMPLIFY_M)
        point_ids: list[int] = []
        for gx, gy in simple:
            point_ids.append(next_point)
            lines.append(f"Point({next_point}) = {{{gx - center[0]:.6f}, {gy - center[1]:.6f}, 0, {lc}}};")
            next_point += 1
        curve_ids: list[int] = []
        for idx, p0 in enumerate(point_ids):
            p1 = point_ids[(idx + 1) % len(point_ids)]
            curve_ids.append(next_curve)
            lines.append(f"Line({next_curve}) = {{{p0}, {p1}}};")
            next_curve += 1
        loop_id = next_loop
        next_loop += 1
        lines.append(f"Curve Loop({loop_id}) = {{{', '.join(str(c) for c in curve_ids)}}};")
        return loop_id

    ring_loops = [add_ring(ring, LAND_MESH_SIZE_M) for ring in shoreline]

    depths: list[int] = []
    for idx, ring in enumerate(shoreline):
        probe = ring[0]
        depth = sum(1 for j, other in enumerate(shoreline) if j != idx and point_in_polygon(probe[0], probe[1], other))
        depths.append(depth)

    island_loops = [loop for loop, depth in zip(ring_loops, depths) if depth % 2 == 0]
    water_surface = 1
    lines.append(f"Plane Surface({water_surface}) = {{1, {', '.join(str(loop) for loop in island_loops)}}};")

    surfaces = [water_surface]
    next_surface = 2
    for idx, (loop, depth) in enumerate(zip(ring_loops, depths)):
        if depth % 2 != 0:
            continue
        hole_loops = []
        for hidx, (hole_loop, hole_depth) in enumerate(zip(ring_loops, depths)):
            if hole_depth == depth + 1 and point_in_polygon(shoreline[hidx][0][0], shoreline[hidx][0][1], shoreline[idx]):
                hole_loops.append(hole_loop)
        surface_loops = [loop] + hole_loops
        lines.append(f"Plane Surface({next_surface}) = {{{', '.join(str(row) for row in surface_loops)}}};")
        surfaces.append(next_surface)
        next_surface += 1

    lines.append(f"Recombine Surface {{{', '.join(str(s) for s in surfaces)}}};")
    for zone in (1, 2, 3):
        lines.append(f"Physical Curve({zone}) = {{{', '.join(str(c) for c in side_curves[zone])}}};")
    lines.append(f"Physical Surface(10) = {{{', '.join(str(s) for s in surfaces)}}};")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_gmsh_geo_uniform(path: Path, radius: float) -> None:
    geo = f"""
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthMin = {MESH_SIZE_M};
Mesh.CharacteristicLengthMax = {MESH_SIZE_M};
Point(1) = {{0, 0, 0, {MESH_SIZE_M}}};
Point(2) = {{{radius}, 0, 0, {MESH_SIZE_M}}};
Point(3) = {{0, {radius}, 0, {MESH_SIZE_M}}};
Point(4) = {{{-radius}, 0, 0, {MESH_SIZE_M}}};
Point(5) = {{0, {-radius}, 0, {MESH_SIZE_M}}};
Circle(1) = {{2, 1, 3}};
Circle(2) = {{3, 1, 4}};
Circle(3) = {{4, 1, 5}};
Circle(4) = {{5, 1, 2}};
Curve Loop(1) = {{1, 2, 3, 4}};
Plane Surface(1) = {{1}};
Recombine Surface {{1}};
Physical Curve(1) = {{2}};
Physical Curve(2) = {{1, 3, 4}};
Physical Surface(10) = {{1}};
"""
    path.write_text(geo.strip() + "\n", encoding="utf-8")


def parse_msh(path: Path) -> tuple[dict[int, tuple[float, float]], list[tuple[int, list[int]]], list[tuple[int, int, int]]]:
    text = path.read_text(encoding="utf-8").splitlines()
    nodes: dict[int, tuple[float, float]] = {}
    quads: list[tuple[int, list[int]]] = []
    lines: list[tuple[int, int, int]] = []
    i = 0
    while i < len(text):
        if text[i] == "$Nodes":
            n = int(text[i + 1])
            for row in text[i + 2 : i + 2 + n]:
                parts = row.split()
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]))
            i += n + 3
            continue
        if text[i] == "$Elements":
            n = int(text[i + 1])
            for row in text[i + 2 : i + 2 + n]:
                p = [int(v) for v in row.split()]
                elm_type, ntags = p[1], p[2]
                tags = p[3 : 3 + ntags]
                conn = p[3 + ntags :]
                physical = tags[0] if tags else 0
                if elm_type == 1:
                    lines.append((physical, conn[0], conn[1]))
                elif elm_type == 3:
                    quads.append((physical, conn))
                elif elm_type in (2,):
                    raise SystemExit("Gmsh produced triangles; expected a recombined quad disk.")
            i += n + 3
            continue
        i += 1
    if not quads:
        raise SystemExit("Gmsh produced no quadrilateral cells.")
    return nodes, quads, lines


def polygon_area(conn: list[int], nodes: dict[int, tuple[float, float]]) -> float:
    pts = [nodes[n] for n in conn]
    return 0.5 * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]))


def write_zone(handle, element: int, face: int, label: int) -> None:
    handle.write(struct.pack("<9i", element, face, 0, label, 0, 0, 0, 0, 7))


def write_nmsh(
    path: Path,
    nodes: dict[int, tuple[float, float]],
    quads: list[tuple[int, list[int]]],
    lines: list[tuple[int, int, int]],
    samples: np.ndarray,
    water_level: float,
    shoreline: list[list[tuple[float, float]]],
    center: tuple[float, float],
    radius: float,
) -> dict:
    base_nodes = sorted(nodes)
    node_pos = {node: idx for idx, node in enumerate(base_nodes)}
    eta = np.linspace(0.0, 1.0, NZ + 1) ** VERTICAL_STRETCH
    bottom = {
        node: terrain_height(x, y, samples, water_level, shoreline, center, radius)
        for node, (x, y) in nodes.items()
    }
    top_z = min(bottom.values()) + DOMAIN_HEIGHT_ABOVE_LOWEST_M

    fixed_quads: list[list[int]] = []
    for _, conn in quads:
        c = list(conn)
        if polygon_area(c, nodes) < 0:
            c.reverse()
        fixed_quads.append(c)

    edge_face: dict[tuple[int, int], tuple[int, int]] = {}
    for qid, conn in enumerate(fixed_quads):
        face_edges = {
            tuple(sorted((conn[0], conn[1]))): 3,
            tuple(sorted((conn[1], conn[2]))): 2,
            tuple(sorted((conn[2], conn[3]))): 4,
            tuple(sorted((conn[3], conn[0]))): 1,
        }
        for edge, face in face_edges.items():
            edge_face[edge] = (qid, face)

    side_zones: list[tuple[int, int, int]] = []
    for label, a, b in lines:
        edge = tuple(sorted((a, b)))
        if edge not in edge_face:
            continue
        qid, face = edge_face[edge]
        zone_label = label if label in (1, 2, 3) else 2
        for ez in range(NZ):
            side_zones.append((1 + qid + len(fixed_quads) * ez, face, zone_label))

    def vertex_global_id(node: int, k: int) -> int:
        return 1 + node_pos[node] + len(base_nodes) * k

    def vertex(node: int, k: int) -> tuple[int, float, float, float]:
        x, y = nodes[node]
        z0 = bottom[node]
        z = z0 + eta[k] * (top_z - z0)
        return vertex_global_id(node, k), x, y, float(z)

    n_elements = len(fixed_quads) * NZ
    zones = list(side_zones)
    for qid in range(len(fixed_quads)):
        zones.append((1 + qid, 5, 5))
        zones.append((1 + qid + len(fixed_quads) * (NZ - 1), 6, 6))

    with path.open("wb") as handle:
        handle.write(struct.pack("<2i", n_elements, 3))
        for ez in range(NZ):
            for qid, conn in enumerate(fixed_quads):
                handle.write(struct.pack("<i", 1 + qid + len(fixed_quads) * ez))
                for idx, x, y, z in (
                    vertex(conn[0], ez),
                    vertex(conn[1], ez),
                    vertex(conn[2], ez),
                    vertex(conn[3], ez),
                    vertex(conn[0], ez + 1),
                    vertex(conn[1], ez + 1),
                    vertex(conn[2], ez + 1),
                    vertex(conn[3], ez + 1),
                ):
                    handle.write(struct.pack("<i3d", idx, x, y, z))
        handle.write(struct.pack("<i", len(zones)))
        for element, face, label in zones:
            write_zone(handle, element, face, label)
        handle.write(struct.pack("<i", 0))

    return {
        "elements": n_elements,
        "base_quads": len(fixed_quads),
        "base_nodes": len(base_nodes),
        "vertices": len(base_nodes) * (NZ + 1),
        "terrain_min": float(min(bottom.values())),
        "terrain_max": float(max(bottom.values())),
        "water_level": water_level,
        "top_z": float(top_z),
        "side_zone_faces": len(side_zones),
    }


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>text{font-family:Arial,sans-serif}.small{font-size:14px}.title{font-size:22px;font-weight:700}</style>",
        '<rect width="100%" height="100%" fill="white"/>',
    ]


def polyline(points: list[tuple[float, float]], project, **attrs) -> str:
    attr = " ".join(f'{k.replace("_", "-")}="{escape(str(v))}"' for k, v in attrs.items())
    pairs = [project(x, y) for x, y in points]
    return f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pairs)}" {attr}/>'


def color_height(z: float, water: float, zmax: float) -> str:
    if z <= water + 0.05:
        return "#9dcfe8"
    t = max(0.0, min(1.0, (z - water) / max(zmax - water, 1.0)))
    stops = [
        (0.00, (190, 226, 170)),
        (0.35, (118, 188, 112)),
        (0.62, (238, 214, 118)),
        (0.82, (190, 135, 75)),
        (1.00, (116, 82, 56)),
    ]
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            a = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
            r = int(c0[0] + a * (c1[0] - c0[0]))
            g = int(c0[1] + a * (c1[1] - c0[1]))
            b = int(c0[2] + a * (c1[2] - c0[2]))
            break
    else:
        r, g, b = stops[-1][1]
    return f"#{r:02x}{g:02x}{b:02x}"


def color_height_rgb(z: float, water: float, zmax: float) -> tuple[float, float, float, float]:
    color = color_height(z, water, zmax)
    return (
        int(color[1:3], 16) / 255.0,
        int(color[3:5], 16) / 255.0,
        int(color[5:7], 16) / 255.0,
        1.0,
    )


def render_terrain_only(
    path: Path,
    shoreline: list[list[tuple[float, float]]],
    contours: Path,
    center: tuple[float, float],
    radius: float,
    samples: np.ndarray,
    water_level: float,
    annotate: bool = True,
) -> None:
    width, height = 1250, 1050
    margin = 90
    scale = min((width - 2 * margin) / (2 * radius), (height - 2 * margin) / (2 * radius))

    def project_global(x: float, y: float) -> tuple[float, float]:
        return width / 2 + (x - center[0]) * scale, height / 2 - (y - center[1]) * scale

    lines = svg_header(width, height)
    if annotate:
        lines.append('<text x="60" y="45" class="title">Södermalm terrain height with contours</text>')
        lines.append(f'<text x="60" y="70" class="small">Color: terrain elevation; contours from source data; water/cylinder edge at z = {water_level:.1f} m</text>')

    n = 130
    zmax = float(np.max(samples[:, 2]))
    cell = 2 * radius / n
    for iy in range(n):
        for ix in range(n):
            x = -radius + (ix + 0.5) * cell
            y = -radius + (iy + 0.5) * cell
            if x * x + y * y > radius * radius:
                continue
            z = terrain_height(x, y, samples, water_level, shoreline, center, radius)
            px, py = project_global(x + center[0] - cell / 2, y + center[1] + cell / 2)
            lines.append(
                f'<rect x="{px:.1f}" y="{py:.1f}" width="{cell * scale + 0.7:.1f}" height="{cell * scale + 0.7:.1f}" '
                f'fill="{color_height(z, water_level, zmax)}" stroke="none"/>'
            )

    for feature in load_geojson(contours)["features"]:
        z = feature["properties"].get("HOJD")
        if z is None:
            continue
        zf = float(z)
        if zf < water_level + 0.5:
            continue
        color = "#2c241f" if int(round(zf)) % 10 == 0 else "#4f4540"
        width_stroke = 1.1 if int(round(zf)) % 10 == 0 else 0.55
        pts = line_points(feature["geometry"])
        if len(pts) >= 2:
            lines.append(polyline(pts, project_global, fill="none", stroke=color, stroke_width=width_stroke, opacity=0.68))

    cx, cy = project_global(center[0], center[1])
    cr = radius * scale
    lines.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{cr:.1f}" fill="none" stroke="#1464a0" stroke-width="2.5"/>')
    for ring in shoreline:
        lines.append(polyline(ring, project_global, fill="none", stroke="#111111", stroke_width=2.2))

    if annotate:
        legend_x, legend_y = 980, 120
        legend_h, legend_w = 240, 24
        for i in range(80):
            t = i / 79
            z = water_level + t * (zmax - water_level)
            y = legend_y + legend_h * (1 - t)
            lines.append(f'<rect x="{legend_x}" y="{y:.1f}" width="{legend_w}" height="{legend_h / 79 + 1:.1f}" fill="{color_height(z, water_level, zmax)}"/>')
        lines.append(f'<rect x="{legend_x}" y="{legend_y}" width="{legend_w}" height="{legend_h}" fill="none" stroke="#333" stroke-width="1"/>')
        lines.append(f'<text x="{legend_x + 34}" y="{legend_y + 6}" class="small">{zmax:.0f} m</text>')
        lines.append(f'<text x="{legend_x + 34}" y="{legend_y + legend_h}" class="small">{water_level:.0f} m</text>')
    lines.append("</svg>")
    path.with_suffix(".svg").write_text("\n".join(lines), encoding="utf-8")
    run([RSVG, str(path.with_suffix(".svg")), "-o", str(path)])


def render_overview(
    path: Path,
    shoreline: list[list[tuple[float, float]]],
    buildings: list[list[tuple[float, float]]],
    center: tuple[float, float],
    radius: float,
    samples: np.ndarray,
    water_level: float,
) -> None:
    width, height = 1250, 1050
    margin = 90
    scale = min((width - 2 * margin) / (2 * radius), (height - 2 * margin) / (2 * radius))

    def project(x: float, y: float) -> tuple[float, float]:
        return width / 2 + (x - center[0]) * scale, height / 2 - (y - center[1]) * scale

    lines = svg_header(width, height)
    lines.append('<text x="60" y="45" class="title">Södermalm terrain blended to water-level cylinder</text>')
    lines.append(f'<text x="60" y="70" class="small">Water level inferred from local contours: z = {water_level:.1f} m; red arc is the {INFLOW_ARC_WIDTH_DEG:.0f} deg inlet centered on wind from {INFLOW_FROM_DEG:.0f} deg</text>')

    n = 95
    zmax = float(np.max(samples[:, 2]))
    cell = 2 * radius / n
    for iy in range(n):
        for ix in range(n):
            x = -radius + (ix + 0.5) * cell
            y = -radius + (iy + 0.5) * cell
            if x * x + y * y > radius * radius:
                continue
            z = terrain_height(x, y, samples, water_level, shoreline, center, radius)
            px, py = project(x + center[0] - cell / 2, y + center[1] + cell / 2)
            lines.append(
                f'<rect x="{px:.1f}" y="{py:.1f}" width="{cell * scale + 0.8:.1f}" height="{cell * scale + 0.8:.1f}" '
                f'fill="{color_height(z, water_level, zmax)}" stroke="none"/>'
            )

    cx, cy = project(center[0], center[1])
    cr = radius * scale
    lines.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{cr:.1f}" fill="none" stroke="#1464a0" stroke-width="3"/>')
    for arc in inlet_arc_polylines(radius):
        pts = [(center[0] + x, center[1] + y) for x, y in arc]
        lines.append(polyline(pts, project, fill="none", stroke="#d7191c", stroke_width=8))
    for ring in buildings[:6000]:
        lines.append(polyline(ring + [ring[0]], project, fill="none", stroke="#555555", stroke_width=0.55, opacity=0.45))
    for ring in shoreline:
        lines.append(polyline(ring, project, fill="none", stroke="#222222", stroke_width=2.5))
    lines.append("</svg>")
    path.with_suffix(".svg").write_text("\n".join(lines), encoding="utf-8")
    run([RSVG, str(path.with_suffix(".svg")), "-o", str(path)])


def render_mesh(
    path: Path,
    nodes: dict[int, tuple[float, float]],
    quads: list[tuple[int, list[int]]],
    radius: float,
    shoreline: list[list[tuple[float, float]]] | None = None,
    center: tuple[float, float] | None = None,
) -> None:
    width, height = 1150, 1050
    margin = 90
    scale = min((width - 2 * margin) / (2 * radius), (height - 2 * margin) / (2 * radius))

    def project(x: float, y: float) -> tuple[float, float]:
        return width / 2 + x * scale, height / 2 - y * scale

    lines = svg_header(width, height)
    if path.stem.endswith("_notext"):
        include_text = False
    else:
        include_text = True
    if include_text:
        lines.append('<text x="60" y="45" class="title">Gmsh quad-disk cylinder mesh footprint</text>')
        lines.append(f'<text x="60" y="70" class="small">{len(quads)} base quads x {NZ} vertical layers; red arc is the {INFLOW_ARC_WIDTH_DEG:.0f} deg inlet zone</text>')
    for _, conn in quads:
        pts = [nodes[n] for n in conn]
        pts.append(pts[0])
        lines.append(polyline(pts, project, fill="none", stroke="#bdbdbd", stroke_width=0.65))
    outer = [(radius * math.cos(t), radius * math.sin(t)) for t in np.linspace(0, 2 * math.pi, 361)]
    lines.append(polyline(outer, project, fill="none", stroke="#1464a0", stroke_width=3))
    if shoreline is not None and center is not None:
        for ring in shoreline:
            local_ring = [(x - center[0], y - center[1]) for x, y in ring]
            lines.append(polyline(local_ring, project, fill="none", stroke="#111111", stroke_width=2.8, opacity=0.92))
    for arc in inlet_arc_polylines(radius):
        lines.append(polyline(arc, project, fill="none", stroke="#d7191c", stroke_width=8))
    lines.append("</svg>")
    path.with_suffix(".svg").write_text("\n".join(lines), encoding="utf-8")
    run([RSVG, str(path.with_suffix(".svg")), "-o", str(path)])


def render_terrain_mesh_3d(
    path: Path,
    nodes: dict[int, tuple[float, float]],
    quads: list[tuple[int, list[int]]],
    lines_1d: list[tuple[int, int, int]],
    samples: np.ndarray,
    water_level: float,
    shoreline: list[list[tuple[float, float]]],
    center: tuple[float, float],
    radius: float,
) -> None:
    width, height = 1800, 1300
    eta = np.linspace(0.0, 1.0, NZ + 1) ** VERTICAL_STRETCH
    bottom = {
        node: terrain_height(x, y, samples, water_level, shoreline, center, radius)
        for node, (x, y) in nodes.items()
    }
    top_z = min(bottom.values()) + DOMAIN_HEIGHT_ABOVE_LOWEST_M

    def z_at(node: int, k: int) -> float:
        return bottom[node] + eta[k] * (top_z - bottom[node])

    def project_raw(x: float, y: float, z: float) -> tuple[float, float]:
        return 0.95 * x, -0.48 * y - 1.15 * (z - water_level)

    point_cloud: list[tuple[float, float]] = []
    for node, (x, y) in nodes.items():
        point_cloud.append(project_raw(x, y, bottom[node]))
        point_cloud.append(project_raw(x, y, top_z))
    min_x = min(x for x, _ in point_cloud)
    max_x = max(x for x, _ in point_cloud)
    min_y = min(y for _, y in point_cloud)
    max_y = max(y for _, y in point_cloud)
    margin = 58
    scale = min((width - 2 * margin) / (max_x - min_x), (height - 2 * margin) / (max_y - min_y))

    def project(x: float, y: float, z: float) -> tuple[float, float]:
        px, py = project_raw(x, y, z)
        return margin + (px - min_x) * scale, margin + (py - min_y) * scale

    def add_line(items: list[dict], p0: tuple[float, float, float], p1: tuple[float, float, float], color: str, width_px: float, opacity: float = 1.0) -> None:
        depth = -0.52 * (p0[0] + p1[0]) * 0.5 + 0.78 * (p0[1] + p1[1]) * 0.5 + 0.08 * (p0[2] + p1[2]) * 0.5
        items.append({"p0": p0, "p1": p1, "color": color, "width": width_px, "opacity": opacity, "depth": depth})

    def draw_line(item: dict) -> str:
        x0, y0 = project(*item["p0"])
        x1, y1 = project(*item["p1"])
        return (
            f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
            f'stroke="{item["color"]}" stroke-width="{item["width"]}" opacity="{item["opacity"]}" stroke-linecap="round"/>'
        )

    surface_edges: set[tuple[int, int]] = set()
    cross_edges: set[tuple[int, int]] = set()
    for _, conn in quads:
        edges = [(conn[0], conn[1]), (conn[1], conn[2]), (conn[2], conn[3]), (conn[3], conn[0])]
        for a, b in edges:
            surface_edges.add(tuple(sorted((a, b))))
        ys = [nodes[n][1] for n in conn]
        if min(ys) <= 0.0 <= max(ys):
            for a, b in edges:
                cross_edges.add(tuple(sorted((a, b))))

    boundary_edges = {tuple(sorted((a, b))): label for label, a, b in lines_1d}
    items: list[dict] = []

    # Pale terrain-following bottom mesh over the full disk.
    for a, b in sorted(surface_edges):
        x0, y0 = nodes[a]
        x1, y1 = nodes[b]
        add_line(items, (x0, y0, bottom[a] + 0.4), (x1, y1, bottom[b] + 0.4), "#4f5a55", 0.55, 0.34)

    # The top cap is drawn lightly, so the extrusion height is readable but not dominant.
    for a, b in sorted(surface_edges):
        x0, y0 = nodes[a]
        x1, y1 = nodes[b]
        add_line(items, (x0, y0, top_z), (x1, y1, top_z), "#6f8ea2", 0.45, 0.18)

    # Inlet side only: avoid drawing the full cylinder boundary so the mesh
    # itself stays readable.
    for (a, b), label in boundary_edges.items():
        if label != 1:
            continue
        x0, y0 = nodes[a]
        x1, y1 = nodes[b]
        for k in range(NZ + 1):
            add_line(items, (x0, y0, z_at(a, k)), (x1, y1, z_at(b, k)), "#d7191c", 1.55, 0.80)
        for node in (a, b):
            x, y = nodes[node]
            add_line(items, (x, y, bottom[node]), (x, y, top_z), "#d7191c", 0.62, 0.36)

    # Central cut strip: this is the important part for seeing the hex layers.
    for a, b in sorted(cross_edges):
        x0, y0 = nodes[a]
        x1, y1 = nodes[b]
        for k in range(NZ + 1):
            add_line(items, (x0, y0, z_at(a, k)), (x1, y1, z_at(b, k)), "#202020", 1.0 if k in (0, NZ) else 0.62, 0.78)
        for node in (a, b):
            x, y = nodes[node]
            add_line(items, (x, y, bottom[node]), (x, y, top_z), "#202020", 0.78, 0.62)

    lines_svg = svg_header(width, height)
    water_disk = [(radius * math.cos(t), radius * math.sin(t), water_level - 0.3) for t in np.linspace(0, 2 * math.pi, 361)]
    water_pts = [project(x, y, z) for x, y, z in water_disk]
    lines_svg.append(
        f'<polygon points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in water_pts)}" '
        'fill="#a9d8ed" stroke="none" opacity="0.78"/>'
    )
    for item in sorted(items, key=lambda row: row["depth"]):
        lines_svg.append(draw_line(item))
    for arc in inlet_arc_polylines(radius):
        lines_svg.append(polyline(arc, lambda x, y: project(x, y, water_level + 1.0), fill="none", stroke="#d7191c", stroke_width=8.5, stroke_linecap="round"))
    lines_svg.append("</svg>")
    path.with_suffix(".svg").write_text("\n".join(lines_svg), encoding="utf-8")
    run([RSVG, str(path.with_suffix(".svg")), "-o", str(path)])


def ring_centroid(ring: list[tuple[float, float]]) -> tuple[float, float]:
    if not ring:
        return (0.0, 0.0)
    pts = ring[:-1] if ring[0] == ring[-1] else ring
    area = 0.0
    cx = 0.0
    cy = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]):
        cross = x0 * y1 - x1 * y0
        area += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    if abs(area) < 1e-9:
        return (sum(x for x, _ in pts) / len(pts), sum(y for _, y in pts) / len(pts))
    return (cx / (3.0 * area), cy / (3.0 * area))


def clean_building_ring(ring: list[tuple[float, float]]) -> list[tuple[float, float]]:
    cleaned: list[tuple[float, float]] = []
    for point in ring:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def building_ring_area(ring: list[tuple[float, float]]) -> float:
    return 0.5 * sum(
        x0 * y1 - x1 * y0
        for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1])
    )


def building_ring_perimeter(ring: list[tuple[float, float]]) -> float:
    return sum(math.hypot(x1 - x0, y1 - y0) for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]))


def building_ring_min_edge(ring: list[tuple[float, float]]) -> float:
    if len(ring) < 2:
        return 0.0
    return min(math.hypot(x1 - x0, y1 - y0) for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]))


def point_in_triangle(
    point: tuple[float, float],
    a: tuple[float, float],
    b: tuple[float, float],
    c: tuple[float, float],
) -> bool:
    px, py = point
    ax, ay = a
    bx, by = b
    cx, cy = c
    v0x, v0y = cx - ax, cy - ay
    v1x, v1y = bx - ax, by - ay
    v2x, v2y = px - ax, py - ay
    denom = v0x * v1y - v1x * v0y
    if abs(denom) < 1e-12:
        return False
    u = (v2x * v1y - v1x * v2y) / denom
    v = (v0x * v2y - v2x * v0y) / denom
    return u >= -1e-12 and v >= -1e-12 and u + v <= 1.0 + 1e-12


def triangulate_building_ring(ring: list[tuple[float, float]]) -> tuple[list[tuple[int, int, int]], bool]:
    ring = clean_building_ring(ring)
    if len(ring) < 3:
        return [], True
    if building_ring_area(ring) < 0.0:
        ring.reverse()

    remaining = list(range(len(ring)))
    triangles: list[tuple[int, int, int]] = []
    guard = 0
    while len(remaining) > 3 and guard < len(ring) * len(ring):
        guard += 1
        clipped = False
        for pos, curr in enumerate(remaining):
            prev_i = remaining[pos - 1]
            next_i = remaining[(pos + 1) % len(remaining)]
            a = ring[prev_i]
            b = ring[curr]
            c = ring[next_i]
            cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            if cross <= 1e-10:
                continue
            if any(
                idx not in (prev_i, curr, next_i) and point_in_triangle(ring[idx], a, b, c)
                for idx in remaining
            ):
                continue
            triangles.append((prev_i, curr, next_i))
            del remaining[pos]
            clipped = True
            break
        if not clipped:
            return [(0, i, i + 1) for i in range(1, len(ring) - 1)], False
    if len(remaining) == 3:
        triangles.append((remaining[0], remaining[1], remaining[2]))
    return triangles, True


def building_base_and_top(properties: dict, terrain_base: float, water_level: float) -> tuple[float, float]:
    mark = properties.get("MARK_Z")
    roof = properties.get("TAK_Z")
    by_height = properties.get("BYGG_H")
    base = max(water_level, terrain_base + BUILDING_BASE_CLEARANCE_M)

    try:
        h = float(roof) - float(mark)
        if BUILDING_MIN_HEIGHT_M <= h <= BUILDING_MAX_HEIGHT_M:
            return base, base + h
    except (TypeError, ValueError):
        pass

    try:
        h = float(by_height)
        if BUILDING_MIN_HEIGHT_M <= h <= BUILDING_MAX_HEIGHT_M:
            return base, base + h
    except (TypeError, ValueError):
        pass
    return base, base + 12.0


def building_sits_on_water_level_terrain(terrain_base: float, water_level: float) -> bool:
    return terrain_base <= water_level + BUILDING_MIN_TERRAIN_ABOVE_WATER_M


def facet_normal(a: tuple[float, float, float], b: tuple[float, float, float], c: tuple[float, float, float]) -> tuple[float, float, float]:
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    length = math.sqrt(nx * nx + ny * ny + nz * nz)
    if length == 0.0:
        return (0.0, 0.0, 0.0)
    return (nx / length, ny / length, nz / length)


def write_sodermalm_building_stl(
    path: Path,
    buildings_geojson: Path,
    samples: np.ndarray,
    water_level: float,
    shoreline: list[list[tuple[float, float]]],
    center: tuple[float, float],
    radius: float,
) -> dict:
    triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]] = []
    skipped = 0
    skipped_small = 0
    skipped_skinny = 0
    fallback_triangulations = 0
    parts = 0

    for feature in load_geojson(buildings_geojson)["features"]:
        props = feature.get("properties", {})
        for ring in polygon_rings(feature["geometry"]):
            pts_global = clean_building_ring(ring)
            if len(pts_global) < 3:
                skipped += 1
                continue
            if BUILDING_SIMPLIFY_M > 0.0:
                pts_global = simplify_ring(pts_global, BUILDING_SIMPLIFY_M)
                pts_global = clean_building_ring(pts_global)
                if len(pts_global) < 3:
                    skipped += 1
                    continue
            footprint_area = abs(building_ring_area(pts_global))
            if footprint_area < BUILDING_MIN_FOOTPRINT_AREA_M2:
                skipped += 1
                skipped_small += 1
                continue
            if building_ring_min_edge(pts_global) < BUILDING_MIN_EDGE_M:
                skipped += 1
                skipped_skinny += 1
                continue
            cx, cy = ring_centroid(pts_global)
            lx, ly = cx - center[0], cy - center[1]
            if lx * lx + ly * ly > radius * radius:
                skipped += 1
                continue
            terrain_base = terrain_height(lx, ly, samples, water_level, shoreline, center, radius)
            if building_sits_on_water_level_terrain(terrain_base, water_level):
                skipped += 1
                continue
            base_z, top_z = building_base_and_top(props, terrain_base, water_level)
            if top_z <= base_z + 0.5:
                skipped += 1
                continue
            pts_local = [(x - center[0], y - center[1]) for x, y in pts_global]
            if building_ring_area(pts_local) < 0.0:
                pts_local.reverse()
            roof_tris, clean = triangulate_building_ring(pts_local)
            if not clean:
                fallback_triangulations += 1

            bottom = [(x, y, base_z) for x, y in pts_local]
            roof = [(x, y, top_z) for x, y in pts_local]
            for ia, ib, ic in roof_tris:
                triangles.append((roof[ia], roof[ib], roof[ic]))
                triangles.append((bottom[ic], bottom[ib], bottom[ia]))
            for i in range(len(pts_local)):
                j = (i + 1) % len(pts_local)
                triangles.append((bottom[i], bottom[j], roof[j]))
                triangles.append((bottom[i], roof[j], roof[i]))
            parts += 1

    with path.open("wb") as handle:
        header = f"Sodermalm sealed building STL: {path.name}".encode("ascii", errors="ignore")[:80]
        handle.write(header.ljust(80, b" "))
        handle.write(struct.pack("<I", len(triangles)))
        for a, b, c in triangles:
            normal = facet_normal(a, b, c)
            handle.write(struct.pack("<12fH", *normal, *a, *b, *c, 0))

    return {
        "path": str(path),
        "facets": len(triangles),
        "building_parts": parts,
        "skipped_parts": skipped,
        "skipped_small_parts": skipped_small,
        "skipped_skinny_parts": skipped_skinny,
        "fallback_triangulations": fallback_triangulations,
        "min_footprint_area_m2": BUILDING_MIN_FOOTPRINT_AREA_M2,
        "min_edge_m": BUILDING_MIN_EDGE_M,
        "simplify_m": BUILDING_SIMPLIFY_M,
        "base_clearance_m": BUILDING_BASE_CLEARANCE_M,
        "height_limits_m": [BUILDING_MIN_HEIGHT_M, BUILDING_MAX_HEIGHT_M],
        "description": "Sealed LOD1 building solids in the local cylinder coordinate system for Brinkman boundary_mesh.",
    }


def shade_color(color: tuple[float, float, float, float], shade: float) -> str:
    shade = max(0.50, min(1.18, shade))
    rgb = [max(0, min(255, int(round(c * 255 * shade)))) for c in color[:3]]
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def render_3d_scene(
    path: Path,
    shoreline: list[list[tuple[float, float]]],
    buildings_geojson: Path,
    center: tuple[float, float],
    radius: float,
    samples: np.ndarray,
    water_level: float,
) -> None:
    width, height = 1800, 1300
    polygons: list[dict] = []
    wall_polygons: list[dict] = []
    wall_edges: list[tuple[tuple[float, float, float], tuple[float, float, float]]] = []
    roof_polygons: list[dict] = []

    def rgb_to_hex(color: tuple[float, float, float, float]) -> str:
        return "#" + "".join(f"{max(0, min(255, int(round(c * 255)))):02x}" for c in color[:3])

    vertical_exaggeration = 1.15

    def project_raw(x: float, y: float, z: float) -> tuple[float, float]:
        # Keep the map close to north-up/east-right while lifting z enough to
        # reveal terrain/buildings without turning it into an isometric rotation.
        return 0.95 * x, -0.48 * y - vertical_exaggeration * (z - water_level)

    def add_poly(
        points: list[tuple[float, float, float]],
        fill: str,
        stroke: str = "none",
        stroke_width: float = 0.0,
        opacity: float = 1.0,
        depth_bias: float = 0.0,
        wall: bool = False,
        roof: bool = False,
    ) -> None:
        if len(points) < 3:
            return
        cx = sum(p[0] for p in points) / len(points)
        cy = sum(p[1] for p in points) / len(points)
        cz = sum(p[2] for p in points) / len(points)
        # Back-to-front painter sorting for the chosen oblique view.
        depth = -0.55 * cx + 0.80 * cy + 0.08 * cz + depth_bias
        item = {"points": points, "fill": fill, "stroke": stroke, "stroke_width": stroke_width, "opacity": opacity, "depth": depth}
        if roof:
            roof_polygons.append(item)
        elif wall:
            wall_polygons.append(item)
        else:
            polygons.append(item)

    n = 104 * 4
    xs = np.linspace(-radius, radius, n + 1)
    ys = np.linspace(-radius, radius, n + 1)
    zmax = float(np.max(samples[:, 2]))
    zgrid = np.full((n + 1, n + 1), water_level, dtype=float)
    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            rr = math.hypot(float(x), float(y))
            sx, sy = (float(x), float(y)) if rr <= radius else (float(x) * radius / rr, float(y) * radius / rr)
            zgrid[iy, ix] = terrain_height(sx, sy, samples, water_level, shoreline, center, radius)
    for iy in range(n):
        for ix in range(n):
            x0, x1 = float(xs[ix]), float(xs[ix + 1])
            y0, y1 = float(ys[iy]), float(ys[iy + 1])
            xc = 0.5 * (x0 + x1)
            yc = 0.5 * (y0 + y1)
            if xc * xc + yc * yc > radius * radius:
                continue
            corners: list[tuple[float, float, float]] = []
            for (x, y), z in zip(
                ((x0, y0), (x1, y0), (x1, y1), (x0, y1)),
                (zgrid[iy, ix], zgrid[iy, ix + 1], zgrid[iy + 1, ix + 1], zgrid[iy + 1, ix]),
            ):
                rr = math.hypot(x, y)
                if rr > radius:
                    x *= radius / rr
                    y *= radius / rr
                corners.append((x, y, float(z)))
            zc = sum(p[2] for p in corners) / 4
            dzdx = ((zgrid[iy, ix + 1] + zgrid[iy + 1, ix + 1]) - (zgrid[iy, ix] + zgrid[iy + 1, ix])) / (2 * (x1 - x0))
            dzdy = ((zgrid[iy + 1, ix] + zgrid[iy + 1, ix + 1]) - (zgrid[iy, ix] + zgrid[iy, ix + 1])) / (2 * (y1 - y0))
            normal = np.array([-dzdx * vertical_exaggeration, -dzdy * vertical_exaggeration, 1.0], dtype=float)
            normal /= max(float(np.linalg.norm(normal)), 1e-12)
            light = np.array([-0.35, -0.55, 0.76], dtype=float)
            light /= float(np.linalg.norm(light))
            hillshade = 0.70 + 0.38 * max(0.0, float(np.dot(normal, light)))
            color = shade_color(color_height_rgb(zc, water_level, zmax), hillshade)
            add_poly(corners, color, opacity=1.0)

    for feature in load_geojson(buildings_geojson)["features"]:
        props = feature.get("properties", {})
        for ring in polygon_rings(feature["geometry"]):
            pts_global = ring[:-1] if ring and ring[0] == ring[-1] else ring
            if len(pts_global) < 3:
                continue
            cx, cy = ring_centroid(pts_global)
            lx, ly = cx - center[0], cy - center[1]
            if lx * lx + ly * ly > radius * radius:
                continue
            pts = [(x - center[0], y - center[1]) for x, y in pts_global]
            terrain_base = terrain_height(lx, ly, samples, water_level, shoreline, center, radius)
            if building_sits_on_water_level_terrain(terrain_base, water_level):
                continue
            base, top = building_base_and_top(props, terrain_base, water_level)
            roof_fill = "#d8d7d2"
            roof_points = [(x, y, top) for x, y in pts]
            for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]):
                dx, dy = x1 - x0, y1 - y0
                length = max(math.hypot(dx, dy), 1e-9)
                normal_x, normal_y = dy / length, -dx / length
                shade = 0.52 + 0.30 * max(0.0, normal_x * -0.45 + normal_y * -0.65) + 0.12 * max(0.0, normal_x * 0.65 + normal_y * 0.20)
                shade = max(0.40, min(0.78, shade))
                c = int(round(shade * 255))
                wall_face = [(x0, y0, base), (x1, y1, base), (x1, y1, top), (x0, y0, top)]
                add_poly(wall_face, f"#{c:02x}{c:02x}{max(0, c - 4):02x}", "#444440", 0.48, wall=True)
                wall_edges.append(((x0, y0, base), (x0, y0, top)))
                wall_edges.append(((x1, y1, base), (x1, y1, top)))
            add_poly(roof_points, roof_fill, "#40403c", 0.75, depth_bias=40.0, roof=True)

    projected = [project_raw(x, y, z) for poly in polygons + wall_polygons + roof_polygons for x, y, z in poly["points"]]
    min_x = min(x for x, _ in projected)
    max_x = max(x for x, _ in projected)
    min_y = min(y for _, y in projected)
    max_y = max(y for _, y in projected)
    margin = 55
    scale = min((width - 2 * margin) / (max_x - min_x), (height - 2 * margin) / (max_y - min_y))

    def project(x: float, y: float, z: float) -> tuple[float, float]:
        px, py = project_raw(x, y, z)
        return margin + (px - min_x) * scale, margin + (py - min_y) * scale

    lines = svg_header(width, height)
    for poly in sorted(polygons, key=lambda item: item["depth"]):
        pts = [project(x, y, z) for x, y, z in poly["points"]]
        stroke = poly["stroke"]
        sw = poly["stroke_width"]
        opacity = poly["opacity"]
        lines.append(
            f'<polygon points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pts)}" '
            f'fill="{poly["fill"]}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>'
        )
    for poly in sorted(wall_polygons, key=lambda item: item["depth"]):
        pts = [project(x, y, z) for x, y, z in poly["points"]]
        lines.append(
            f'<polygon points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pts)}" '
            f'fill="{poly["fill"]}" stroke="{poly["stroke"]}" stroke-width="{poly["stroke_width"]}" opacity="{poly["opacity"]}"/>'
        )
    for p0, p1 in wall_edges:
        x0, y0 = project(*p0)
        x1, y1 = project(*p1)
        lines.append(f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" stroke="#3f3f3b" stroke-width="0.35" opacity="0.45"/>')
    for poly in sorted(roof_polygons, key=lambda item: item["depth"]):
        pts = [project(x, y, z) for x, y, z in poly["points"]]
        lines.append(
            f'<polygon points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pts)}" '
            f'fill="{poly["fill"]}" stroke="{poly["stroke"]}" stroke-width="{poly["stroke_width"]}" opacity="{poly["opacity"]}"/>'
        )
    outer = [(radius * math.cos(t), radius * math.sin(t), water_level + 0.05) for t in np.linspace(0, 2 * math.pi, 361)]
    lines.append(polyline([(x, y) for x, y, _ in outer], lambda x, y: project(x, y, water_level + 0.05), fill="none", stroke="#0b5f9f", stroke_width=3.2))
    for arc in inlet_arc_polylines(radius):
        lines.append(polyline(arc, lambda x, y: project(x, y, water_level + 1.0), fill="none", stroke="#d7191c", stroke_width=8.5, stroke_linecap="round"))
    lines.append("</svg>")
    path.with_suffix(".svg").write_text("\n".join(lines), encoding="utf-8")
    run([RSVG, str(path.with_suffix(".svg")), "-o", str(path)])


def render_oblique(path: Path, radius: float, terrain_min: float, top_z: float) -> None:
    width, height = 1200, 850

    def project(x: float, y: float, z: float) -> tuple[float, float]:
        return width / 2 + 0.062 * (x - y), 610 + 0.026 * (x + y) - 1.35 * (z - terrain_min)

    lines = svg_header(width, height)
    lines.append('<text x="60" y="45" class="title">Oblique check: cylindrical side, flat water edge, flat top</text>')
    lines.append(f'<text x="60" y="70" class="small">Blue: cylinder wall; red: {INFLOW_ARC_WIDTH_DEG:.0f} deg inlet; bottom circle is at inferred water level</text>')
    for z, stroke in ((terrain_min, "#1464a0"), (top_z, "#1464a0")):
        pts = [(radius * math.cos(t), radius * math.sin(t)) for t in np.linspace(0, 2 * math.pi, 160)]
        lines.append(polyline(pts, lambda x, y, zz=z: project(x, y, zz), fill="none", stroke=stroke, stroke_width=2))
    for a in np.linspace(0, 360, 17):
        x = radius * math.cos(math.radians(a))
        y = radius * math.sin(math.radians(a))
        p0 = project(x, y, terrain_min)
        p1 = project(x, y, top_z)
        lines.append(f'<line x1="{p0[0]:.1f}" y1="{p0[1]:.1f}" x2="{p1[0]:.1f}" y2="{p1[1]:.1f}" stroke="#9ecae1" stroke-width="1"/>')
    for arc in inlet_arc_polylines(radius):
        lines.append(polyline(arc, lambda x, y: project(x, y, top_z), fill="none", stroke="#d7191c", stroke_width=7))
    lines.append("</svg>")
    path.with_suffix(".svg").write_text("\n".join(lines), encoding="utf-8")
    run([RSVG, str(path.with_suffix(".svg")), "-o", str(path)])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    if not INPUT_GPKG.exists():
        raise SystemExit(f"Missing input GeoPackage: {INPUT_GPKG}")

    run([OSGEO_PYTHON, str(HERE / "fetch_sodermalm_osm_mask.py")])
    land_mask = OUT / "sodermalm_osm_island_epsg3006.geojson"
    land_geometry = load_geojson(land_mask)["features"][0]["geometry"]
    shoreline = polygon_boundary_rings(land_geometry)
    all_shore_points = [p for ring in shoreline for p in ring]
    px = np.array([p[0] for p in all_shore_points])
    py = np.array([p[1] for p in all_shore_points])
    center = (float((px.min() + px.max()) * 0.5), float((py.min() + py.max()) * 0.5))
    radius = float(max(math.hypot(x - center[0], y - center[1]) for x, y in all_shore_points) + BUFFER_M)

    buildings_geojson = OUT / "sodermalm_buildings.geojson"
    contours_geojson = OUT / "sodermalm_cylinder_contours.geojson"
    for path in (buildings_geojson, contours_geojson):
        if path.exists():
            path.unlink()
    run([OGR2OGR, "-f", "GeoJSON", str(buildings_geojson), str(INPUT_GPKG), "buildings", "-clipsrc", str(land_mask)])
    run([
        OGR2OGR, "-f", "GeoJSON", str(contours_geojson), str(INPUT_GPKG), "terrain_contours",
        "-spat", str(center[0] - radius), str(center[1] - radius), str(center[0] + radius), str(center[1] + radius),
    ])

    samples, water_level = terrain_samples(contours_geojson, center)
    geo = OUT / "sodermalm_cylinder_quad_disk.geo"
    msh = OUT / "sodermalm_cylinder_quad_disk.msh"
    write_gmsh_geo(geo, radius, shoreline, center)
    run([GMSH, "-2", str(geo), "-format", "msh2", "-o", str(msh)])
    nodes, quads, lines = parse_msh(msh)

    mesh_path = OUT / "sodermalm_cylinder_coarse.nmsh"
    mesh_info = write_nmsh(mesh_path, nodes, quads, lines, samples, water_level, shoreline, center, radius)
    building_stl_info = write_sodermalm_building_stl(
        OUT / "sodermalm_buildings.stl",
        buildings_geojson,
        samples,
        water_level,
        shoreline,
        center,
        radius,
    )

    building_rings: list[list[tuple[float, float]]] = []
    for feature in load_geojson(buildings_geojson)["features"]:
        for ring in polygon_rings(feature["geometry"]):
            if ring:
                building_rings.append(ring)

    render_overview(FIG / "sodermalm_waterlevel_overview.png", shoreline, building_rings, center, radius, samples, water_level)
    render_terrain_only(FIG / "sodermalm_terrain_height_contours.png", shoreline, contours_geojson, center, radius, samples, water_level)
    render_terrain_only(FIG / "sodermalm_terrain_height_contours_notext.png", shoreline, contours_geojson, center, radius, samples, water_level, annotate=False)
    render_mesh(FIG / "sodermalm_cylinder_mesh_topdown.png", nodes, quads, radius, shoreline, center)
    render_mesh(FIG / "sodermalm_cylinder_mesh_topdown_notext.png", nodes, quads, radius, shoreline, center)
    render_terrain_mesh_3d(FIG / "sodermalm_cylinder_mesh_3d.png", nodes, quads, lines, samples, water_level, shoreline, center, radius)
    render_oblique(FIG / "sodermalm_cylinder_oblique.png", radius, water_level, mesh_info["top_z"])
    render_3d_scene(FIG / "sodermalm_cylinder_3d_buildings.png", shoreline, buildings_geojson, center, radius, samples, water_level)
    checker_output = run_mesh_checker(mesh_path)

    metadata = {
        "center_epsg3006": center,
        "radius_m": radius,
        "buffer_m": BUFFER_M,
        "domain_height_above_lowest_m": DOMAIN_HEIGHT_ABOVE_LOWEST_M,
        "mesh_size_m": MESH_SIZE_M,
        "land_mesh_size_m": LAND_MESH_SIZE_M,
        "water_mesh_size_m": WATER_MESH_SIZE_M,
        "shoreline_simplify_m": SHORELINE_SIMPLIFY_M,
        "shore_blend_m": SHORE_BLEND_M,
        "inflow_from_degrees": INFLOW_FROM_DEG,
        "inflow_arc_width_degrees": INFLOW_ARC_WIDTH_DEG,
        "outflow_arc_width_degrees": OUTFLOW_ARC_WIDTH_DEG,
        "inflow_angle_convention": "meteorological: north is 0/360 degrees, clockwise positive; direction is where wind comes from",
        "zones": {
            "1": "inlet_arc",
            "2": "downstream_outlet_arc",
            "3": "lateral_open_side_arcs",
            "5": "terrain_or_water_bottom",
            "6": "top",
        },
        "mesh": mesh_info,
        "building_stl": building_stl_info,
        "building_footprints": len(building_rings),
        "land_mask": "OpenStreetMap place=island polygon for Södermalm",
        "land_mask_source": str(land_mask),
        "rough_search_hint_wgs84": SODERMALM_WGS84,
        "mesh_checker_summary": checker_output,
    }
    (OUT / "sodermalm_cylinder_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in metadata.items() if k != "mesh_checker_summary"}, indent=2))


if __name__ == "__main__":
    main()
