#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import os
import shutil
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import shapely
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
INPUT_GPKG = PACKAGE_ROOT / "input" / "steep_pilot_layers.gpkg"
OUT_DIR = PACKAGE_ROOT / "generated"
BUILDING_LAYER = "buildings"
TERRAIN_LAYER = "terrain_contours"
PILOT_BBOX = (676000.0, 6577000.0, 678000.0, 6579000.0)

OUTPUT_SLUG = "steep_pilot"
OUTPUT_LABEL = "Steep Pilot DEM-fixed"
GRID_SIZE = 201
CONTOUR_DENSIFY_M = 5.0
TERRAIN_SMOOTHING_M = 20.0
SUSPICIOUS_HEIGHT_M = 80.0
CARVE_TERRAIN_UNDER_BUILDINGS = True
BOX_TOP_Z = 2000.0
STITCHED_EDGE_M = 10.0
TERRAIN_BOTTOM_MARGIN_M = 50.0


@dataclass
class Mesh:
    name: str
    vertices: list[tuple[float, float, float]]
    facets: list[tuple[int, int, int]]


@dataclass
class TerrainSampler:
    xs: np.ndarray
    ys: np.ndarray
    z: np.ndarray

    def sample(self, x: float, y: float) -> float:
        i = int(np.searchsorted(self.xs, x) - 1)
        j = int(np.searchsorted(self.ys, y) - 1)
        i = max(0, min(i, len(self.xs) - 2))
        j = max(0, min(j, len(self.ys) - 2))
        x0, x1 = float(self.xs[i]), float(self.xs[i + 1])
        y0, y1 = float(self.ys[j]), float(self.ys[j + 1])
        tx = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
        ty = 0.0 if y1 == y0 else (y - y0) / (y1 - y0)
        z00 = float(self.z[j, i])
        z10 = float(self.z[j, i + 1])
        z01 = float(self.z[j + 1, i])
        z11 = float(self.z[j + 1, i + 1])
        return (1 - tx) * (1 - ty) * z00 + tx * (1 - ty) * z10 + (1 - tx) * ty * z01 + tx * ty * z11


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_capture(cmd: list[str]) -> str:
    print("+", " ".join(cmd))
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    if result.stderr:
        print(result.stderr, end="")
    return result.stdout


def require_tools(names: Iterable[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise SystemExit(f"Missing required tools: {', '.join(missing)}")


def export_layer(gpkg: Path, layer: str, out_geojson: Path) -> None:
    if out_geojson.exists():
        out_geojson.unlink()
    run(
        [
            "ogr2ogr",
            "-f",
            "GeoJSON",
            str(out_geojson),
            str(gpkg),
            layer,
        ]
    )


def load_geojson(path: Path) -> list[dict]:
    return json.loads(path.read_text(encoding="utf-8"))["features"]


def rings_from_geometry(geometry: dict) -> list[list[list[tuple[float, float]]]]:
    if not geometry:
        return []
    gtype = geometry["type"]
    coords = geometry["coordinates"]
    if gtype == "Polygon":
        return [[[(float(x), float(y)) for x, y, *_ in ring] for ring in coords]]
    if gtype == "MultiPolygon":
        return [[[(float(x), float(y)) for x, y, *_ in ring] for ring in polygon] for polygon in coords]
    return []


def geometry_bounds(geometry: dict) -> tuple[float, float, float, float] | None:
    points: list[tuple[float, float]] = []
    for polygon in rings_from_geometry(geometry):
        for ring in polygon:
            points.extend(ring)
    for part in line_parts(geometry):
        points.extend(part)
    if not points:
        return None
    xs = [x for x, _ in points]
    ys = [y for _, y in points]
    return min(xs), min(ys), max(xs), max(ys)


def infer_bbox(features: list[dict], padding_m: float = 0.0) -> tuple[float, float, float, float]:
    bounds = [geometry_bounds(feature.get("geometry")) for feature in features]
    bounds = [bound for bound in bounds if bound is not None]
    if not bounds:
        raise SystemExit("Could not infer --bbox because input layers contain no geometry.")
    xmin = min(bound[0] for bound in bounds) - padding_m
    ymin = min(bound[1] for bound in bounds) - padding_m
    xmax = max(bound[2] for bound in bounds) + padding_m
    ymax = max(bound[3] for bound in bounds) + padding_m
    if xmax <= xmin or ymax <= ymin:
        raise SystemExit(f"Invalid inferred bbox: {(xmin, ymin, xmax, ymax)}")
    return (xmin, ymin, xmax, ymax)


def line_parts(geometry: dict) -> list[list[tuple[float, float]]]:
    if not geometry:
        return []
    gtype = geometry["type"]
    coords = geometry["coordinates"]
    if gtype == "LineString":
        return [[(float(x), float(y)) for x, y, *_ in coords]]
    if gtype == "MultiLineString":
        return [[(float(x), float(y)) for x, y, *_ in part] for part in coords]
    return []


def clean_ring(ring: list[tuple[float, float]]) -> list[tuple[float, float]]:
    cleaned: list[tuple[float, float]] = []
    for point in ring:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def densify_ring(ring: list[tuple[float, float]], max_segment_m: float = STITCHED_EDGE_M) -> list[tuple[float, float]]:
    ring = clean_ring(ring)
    if len(ring) < 2:
        return ring
    dense: list[tuple[float, float]] = []
    for p0, p1 in zip(ring, ring[1:] + ring[:1]):
        x0, y0 = p0
        x1, y1 = p1
        dist = math.hypot(x1 - x0, y1 - y0)
        steps = max(1, int(math.ceil(dist / max_segment_m)))
        for step in range(steps):
            t = step / steps
            dense.append((x0 + t * (x1 - x0), y0 + t * (y1 - y0)))
    return clean_ring(dense)


def polygon_area(ring: list[tuple[float, float]]) -> float:
    return 0.5 * sum(
        x0 * y1 - x1 * y0
        for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1])
    )


def is_point_in_triangle(
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


def triangulate_ring(ring: list[tuple[float, float]]) -> tuple[list[tuple[int, int, int]], bool]:
    ring = clean_ring(ring)
    if len(ring) < 3:
        return [], True
    if polygon_area(ring) < 0:
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
                idx not in (prev_i, curr, next_i) and is_point_in_triangle(ring[idx], a, b, c)
                for idx in remaining
            ):
                continue
            triangles.append((prev_i, curr, next_i))
            del remaining[pos]
            clipped = True
            break
        if not clipped:
            # A fallback fan is acceptable for this preliminary geometry, and
            # the metadata flags any footprint where ear clipping was not clean.
            return [(0, i, i + 1) for i in range(1, len(ring) - 1)], False

    if len(remaining) == 3:
        triangles.append((remaining[0], remaining[1], remaining[2]))
    return triangles, True


def densified_contour_points(features: list[dict]) -> list[tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    for feature in features:
        z = float(feature["properties"]["HOJD"])
        for part in line_parts(feature["geometry"]):
            for p0, p1 in zip(part, part[1:]):
                x0, y0 = p0
                x1, y1 = p1
                dist = math.hypot(x1 - x0, y1 - y0)
                steps = max(1, int(math.ceil(dist / CONTOUR_DENSIFY_M)))
                for step in range(steps):
                    t = step / steps
                    points.append((x0 + t * (x1 - x0), y0 + t * (y1 - y0), z))
            if part:
                points.append((part[-1][0], part[-1][1], z))
    return points


def write_contour_vrt(points: list[tuple[float, float, float]], out_dir: Path) -> Path:
    csv_path = out_dir / "terrain_contour_points.csv"
    vrt_path = out_dir / "terrain_contour_points.vrt"
    with csv_path.open("w", encoding="utf-8") as f:
        f.write("x,y,z\n")
        for x, y, z in points:
            f.write(f"{x:.6f},{y:.6f},{z:.6f}\n")

    vrt_path.write_text(
        f"""<OGRVRTDataSource>
  <OGRVRTLayer name="terrain_contour_points">
    <SrcDataSource relativeToVRT="1">{csv_path.name}</SrcDataSource>
    <GeometryType>wkbPoint25D</GeometryType>
    <LayerSRS>EPSG:3006</LayerSRS>
    <GeometryField encoding="PointFromColumns" x="x" y="y" z="z"/>
  </OGRVRTLayer>
</OGRVRTDataSource>
""",
        encoding="utf-8",
    )
    return vrt_path


def grid_terrain(vrt_path: Path, out_dir: Path) -> tuple[Path, Path]:
    xmin, ymin, xmax, ymax = PILOT_BBOX
    tif_path = out_dir / "terrain_dem.tif"
    xyz_path = out_dir / "terrain_dem.xyz"
    if tif_path.exists():
        tif_path.unlink()
    run(
        [
            "gdal_grid",
            "--quiet",
            "-zfield",
            "z",
            "-txe",
            str(xmin),
            str(xmax),
            "-tye",
            str(ymin),
            str(ymax),
            "-outsize",
            str(GRID_SIZE),
            str(GRID_SIZE),
            "-a",
            f"invdist:power=2.0:smoothing={TERRAIN_SMOOTHING_M}",
            "-of",
            "GTiff",
            str(vrt_path),
            str(tif_path),
        ]
    )
    if xyz_path.exists():
        xyz_path.unlink()
    run(["gdal_translate", "--quiet", "-of", "XYZ", str(tif_path), str(xyz_path)])
    return tif_path, xyz_path


def load_dem_xyz(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    xs = np.unique(data[:, 0])
    ys = np.unique(data[:, 1])
    z = np.empty((len(ys), len(xs)))
    for x, y, value in data:
        i = int(np.searchsorted(xs, x))
        j = int(np.searchsorted(ys, y))
        z[j, i] = value
    return xs, ys, z


def terrain_mesh(xs: np.ndarray, ys: np.ndarray, z: np.ndarray) -> Mesh:
    xmin, ymin, _, _ = PILOT_BBOX
    vertices: list[tuple[float, float, float]] = []
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            vertices.append((float(x - xmin), float(y - ymin), float(z[j, i])))

    facets: list[tuple[int, int, int]] = []
    nx = len(xs)
    ny = len(ys)
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            b = j * nx + i + 1
            c = (j + 1) * nx + i
            d = (j + 1) * nx + i + 1
            facets.append((a, b, d))
            facets.append((a, d, c))
    return Mesh("terrain", vertices, facets)


def footprint_polygons(features: list[dict]) -> list[Polygon]:
    polygons: list[Polygon] = []
    for feature in features:
        for polygon in rings_from_geometry(feature["geometry"]):
            if not polygon:
                continue
            outer = clean_ring(polygon[0])
            if len(outer) < 3:
                continue
            holes = [clean_ring(ring) for ring in polygon[1:] if len(clean_ring(ring)) >= 3]
            poly = Polygon(outer, holes)
            if not poly.is_empty and poly.area > 0:
                if not poly.is_valid:
                    poly = poly.buffer(0)
                if not poly.is_empty:
                    polygons.append(poly)
    return polygons


def footprint_polygons_densified(features: list[dict]) -> list[Polygon]:
    polygons: list[Polygon] = []
    for feature in features:
        for polygon in rings_from_geometry(feature["geometry"]):
            if not polygon:
                continue
            outer = densify_ring(polygon[0])
            if len(outer) < 3:
                continue
            holes = [densify_ring(ring) for ring in polygon[1:] if len(clean_ring(ring)) >= 3]
            poly = Polygon(outer, holes)
            if not poly.is_empty and poly.area > 0:
                if not poly.is_valid:
                    poly = poly.buffer(0)
                if not poly.is_empty:
                    polygons.append(poly)
    return polygons


def building_union_regions(features: list[dict]) -> tuple[list[Polygon], list[float], dict]:
    source_polygons: list[Polygon] = []
    source_heights: list[float] = []
    skipped = 0
    for feature in features:
        props = feature["properties"]
        mark_z = numeric(props, "MARK_Z")
        top_z = numeric(props, "TAK_Z")
        bygg_h = numeric(props, "BYGG_H")
        if mark_z is None:
            skipped += 1
            continue
        if top_z is None and bygg_h is not None:
            top_z = mark_z + bygg_h
        if top_z is None or top_z <= mark_z:
            skipped += 1
            continue
        height = top_z - mark_z
        for polygon in rings_from_geometry(feature["geometry"]):
            if not polygon:
                continue
            outer = clean_ring(polygon[0])
            if len(outer) < 3:
                continue
            poly = Polygon(outer)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if not poly.is_empty and poly.area > 0:
                source_polygons.append(poly)
                source_heights.append(height)

    if not source_polygons:
        return [], [], {"source_polygons": 0, "regions": 0, "skipped": skipped}

    unioned = unary_union(source_polygons)
    tree = STRtree(source_polygons)
    regions: list[Polygon] = []
    heights: list[float] = []
    for region in iter_polygons(unioned):
        if region.is_empty or region.area <= 0:
            continue
        probe = region.representative_point()
        candidates = [
            int(index)
            for index in tree.query(probe)
            if source_polygons[int(index)].intersects(region)
        ]
        height = max((source_heights[index] for index in candidates), default=3.0)
        outer = densify_ring(list(region.exterior.coords)[:-1])
        holes = [densify_ring(list(interior.coords)[:-1]) for interior in region.interiors]
        dense_region = Polygon(outer, holes)
        if not dense_region.is_valid:
            dense_region = dense_region.buffer(0)
        for dense_part in iter_polygons(dense_region):
            if not dense_part.is_empty and dense_part.area > 0:
                regions.append(dense_part)
                heights.append(height)

    return regions, heights, {
        "source_polygons": len(source_polygons),
        "regions": len(regions),
        "skipped": skipped,
    }


def write_cgal_input(
    path: Path,
    xs: np.ndarray,
    ys: np.ndarray,
    z: np.ndarray,
    regions: list[Polygon],
    heights: list[float],
) -> None:
    xmin, ymin, xmax, ymax = PILOT_BBOX
    with path.open("w", encoding="utf-8") as f:
        f.write(f"bbox {xmin:.12g} {ymin:.12g} {xmax:.12g} {ymax:.12g}\n")
        f.write(f"grid {len(xs)} {len(ys)}\n")
        f.write("xs " + " ".join(f"{float(x):.12g}" for x in xs) + "\n")
        f.write("ys " + " ".join(f"{float(y):.12g}" for y in ys) + "\n")
        f.write("z\n")
        for value in z.ravel():
            f.write(f"{float(value):.12g} ")
        f.write("\n")
        f.write(f"buildings {len(regions)}\n")
        for region, height in zip(regions, heights):
            ring = clean_ring([(float(x), float(y)) for x, y in region.exterior.coords])
            if len(ring) < 3:
                continue
            f.write(f"height {float(height):.12g}\n")
            f.write(f"ring {len(ring)}\n")
            for x, y in ring:
                f.write(f"{x:.12g} {y:.12g}\n")


def filter_interior_building_regions(
    regions: list[Polygon],
    heights: list[float],
    eps_m: float = 1e-6,
    simplify_tolerance_m: float = 0.05,
) -> tuple[list[Polygon], list[float], dict]:
    xmin, ymin, xmax, ymax = PILOT_BBOX
    kept_regions: list[Polygon] = []
    kept_heights: list[float] = []
    skipped_boundary = 0
    for region, height in zip(regions, heights):
        minx, miny, maxx, maxy = region.bounds
        touches_boundary = (
            minx <= xmin + eps_m
            or miny <= ymin + eps_m
            or maxx >= xmax - eps_m
            or maxy >= ymax - eps_m
        )
        if touches_boundary:
            skipped_boundary += 1
            continue
        simplified = region.simplify(simplify_tolerance_m, preserve_topology=True)
        for part in iter_polygons(simplified):
            if not part.is_empty and part.area > 0:
                kept_regions.append(part)
                kept_heights.append(height)
    return kept_regions, kept_heights, {
        "skipped_boundary_clipped_regions": skipped_boundary,
        "kept_regions": len(kept_regions),
        "simplify_tolerance_m": simplify_tolerance_m,
    }


def build_cgal_helper() -> Path:
    source_dir = PACKAGE_ROOT / "src" / "urban_cgal_stl"
    build_dir = PACKAGE_ROOT / ".build" / "urban_cgal_stl"
    build_dir.mkdir(parents=True, exist_ok=True)
    run(["cmake", "-S", str(source_dir), "-B", str(build_dir), "-DCMAKE_BUILD_TYPE=Release"])
    run(["cmake", "--build", str(build_dir), "--config", "Release", "-j"])
    exe = build_dir / "urban_cgal_stl"
    if not exe.exists():
        raise SystemExit(f"CGAL helper build did not produce {exe}")
    return exe


def run_cgal_stl(
    out_dir: Path,
    xs: np.ndarray,
    ys: np.ndarray,
    z: np.ndarray,
    regions: list[Polygon],
    heights: list[float],
    top_z: float,
) -> dict:
    helper = build_cgal_helper()
    cgal_input = out_dir / "urban_cgal_input.txt"
    write_cgal_input(cgal_input, xs, ys, z, regions, heights)
    stdout = run_capture(
        [
            str(helper),
            str(cgal_input),
            str(out_dir / "terrain.stl"),
            str(out_dir / "buildings.stl"),
            str(out_dir / f"{OUTPUT_SLUG}_combined.stl"),
            str(out_dir / f"{OUTPUT_SLUG}_boxed.stl"),
            str(top_z),
        ]
    )
    return json.loads(stdout)


def terrain_domain_with_building_holes(polygons: list[Polygon]) -> Polygon:
    xmin, ymin, xmax, ymax = PILOT_BBOX
    domain = Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])
    if not polygons:
        return domain
    footprints = unary_union(polygons)
    carved = domain.difference(footprints)
    if carved.is_empty:
        raise SystemExit("Building footprints removed the entire terrain domain.")
    return carved


def iter_polygons(geometry) -> Iterable[Polygon]:
    if geometry.geom_type == "Polygon":
        yield geometry
    elif geometry.geom_type == "MultiPolygon":
        yield from geometry.geoms


def stitched_terrain_mesh(domain, sampler: TerrainSampler) -> Mesh:
    xmin, ymin, _, _ = PILOT_BBOX
    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []
    vertex_index: dict[tuple[float, float, float], int] = {}

    def add_vertex(x: float, y: float) -> int:
        z = sampler.sample(x, y)
        key = (round(x - xmin, 6), round(y - ymin, 6), round(z, 6))
        if key not in vertex_index:
            vertex_index[key] = len(vertices)
            vertices.append(key)
        return vertex_index[key]

    triangles = shapely.constrained_delaunay_triangles(domain)
    for tri in triangles.geoms:
        if tri.is_empty or not domain.covers(tri.representative_point()):
            continue
        coords = list(tri.exterior.coords)[:-1]
        if len(coords) != 3:
            continue
        a = add_vertex(float(coords[0][0]), float(coords[0][1]))
        b = add_vertex(float(coords[1][0]), float(coords[1][1]))
        c = add_vertex(float(coords[2][0]), float(coords[2][1]))
        facets.append((a, b, c))
    return Mesh("stitched_terrain", vertices, facets)


def terrain_mesh_carved(xs: np.ndarray, ys: np.ndarray, z: np.ndarray, polygons: list[Polygon]) -> tuple[Mesh, int]:
    xmin, ymin, _, _ = PILOT_BBOX
    vertices: list[tuple[float, float, float]] = []
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            vertices.append((float(x - xmin), float(y - ymin), float(z[j, i])))

    tree = STRtree(polygons)
    facets: list[tuple[int, int, int]] = []
    carved_cells = 0
    nx = len(xs)
    ny = len(ys)
    for j in range(ny - 1):
        for i in range(nx - 1):
            cx = float((xs[i] + xs[i + 1]) * 0.5)
            cy = float((ys[j] + ys[j + 1]) * 0.5)
            point = Point(cx, cy)
            candidate_indexes = tree.query(point)
            if any(polygons[int(index)].contains(point) for index in candidate_indexes):
                carved_cells += 1
                continue
            a = j * nx + i
            b = j * nx + i + 1
            c = (j + 1) * nx + i
            d = (j + 1) * nx + i + 1
            facets.append((a, b, d))
            facets.append((a, d, c))
    return Mesh("terrain_carved", vertices, facets), carved_cells


def terrain_box_mesh(xs: np.ndarray, ys: np.ndarray, z: np.ndarray, top_z: float) -> Mesh:
    xmin, ymin, _, _ = PILOT_BBOX
    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []

    def add_vertex(x: float, y: float, vz: float) -> int:
        vertices.append((float(x - xmin), float(y - ymin), float(vz)))
        return len(vertices) - 1

    def add_quad(a: int, b: int, c: int, d: int) -> None:
        facets.append((a, b, c))
        facets.append((a, c, d))

    nx = len(xs)
    ny = len(ys)
    x0 = float(xs[0])
    x1 = float(xs[-1])
    y0 = float(ys[0])
    y1 = float(ys[-1])

    # Top cap uses the same grid as the terrain so perimeter edges match the
    # segmented side walls exactly.
    top_start = len(vertices)
    for y in ys:
        for x in xs:
            add_vertex(float(x), float(y), top_z)
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = top_start + j * nx + i
            b = top_start + j * nx + i + 1
            c = top_start + (j + 1) * nx + i
            d = top_start + (j + 1) * nx + i + 1
            facets.append((a, d, b))
            facets.append((a, c, d))

    # Side walls follow the terrain along the DEM perimeter.
    for i in range(nx - 1):
        xb0 = float(xs[i])
        xb1 = float(xs[i + 1])
        # South wall.
        b0 = add_vertex(xb0, y0, float(z[0, i]))
        b1 = add_vertex(xb1, y0, float(z[0, i + 1]))
        t1 = add_vertex(xb1, y0, top_z)
        t0 = add_vertex(xb0, y0, top_z)
        add_quad(b0, b1, t1, t0)
        # North wall.
        b0 = add_vertex(xb0, y1, float(z[ny - 1, i]))
        b1 = add_vertex(xb1, y1, float(z[ny - 1, i + 1]))
        t1 = add_vertex(xb1, y1, top_z)
        t0 = add_vertex(xb0, y1, top_z)
        add_quad(b1, b0, t0, t1)

    for j in range(ny - 1):
        yb0 = float(ys[j])
        yb1 = float(ys[j + 1])
        # West wall.
        b0 = add_vertex(x0, yb0, float(z[j, 0]))
        b1 = add_vertex(x0, yb1, float(z[j + 1, 0]))
        t1 = add_vertex(x0, yb1, top_z)
        t0 = add_vertex(x0, yb0, top_z)
        add_quad(b1, b0, t0, t1)
        # East wall.
        b0 = add_vertex(x1, yb0, float(z[j, nx - 1]))
        b1 = add_vertex(x1, yb1, float(z[j + 1, nx - 1]))
        t1 = add_vertex(x1, yb1, top_z)
        t0 = add_vertex(x1, yb0, top_z)
        add_quad(b0, b1, t1, t0)

    return Mesh("box_sides_top", vertices, facets)


def terrain_solid_box_mesh(xs: np.ndarray, ys: np.ndarray, z: np.ndarray, bottom_z: float) -> Mesh:
    xmin, ymin, _, _ = PILOT_BBOX
    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []

    def add_vertex(x: float, y: float, vz: float) -> int:
        vertices.append((float(x - xmin), float(y - ymin), float(vz)))
        return len(vertices) - 1

    def add_quad(a: int, b: int, c: int, d: int) -> None:
        facets.append((a, b, c))
        facets.append((a, c, d))

    nx = len(xs)
    ny = len(ys)

    # Top terrain surface.
    top_start = len(vertices)
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            add_vertex(float(x), float(y), float(z[j, i]))
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = top_start + j * nx + i
            b = top_start + j * nx + i + 1
            c = top_start + (j + 1) * nx + i
            d = top_start + (j + 1) * nx + i + 1
            facets.append((a, b, d))
            facets.append((a, d, c))

    # Flat bottom cap.
    bottom_start = len(vertices)
    for y in ys:
        for x in xs:
            add_vertex(float(x), float(y), bottom_z)
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = bottom_start + j * nx + i
            b = bottom_start + j * nx + i + 1
            c = bottom_start + (j + 1) * nx + i
            d = bottom_start + (j + 1) * nx + i + 1
            facets.append((a, d, b))
            facets.append((a, c, d))

    # Side walls.
    for i in range(nx - 1):
        top0 = top_start + i
        top1 = top_start + i + 1
        bot0 = bottom_start + i
        bot1 = bottom_start + i + 1
        add_quad(bot0, bot1, top1, top0)

        top0 = top_start + (ny - 1) * nx + i
        top1 = top_start + (ny - 1) * nx + i + 1
        bot0 = bottom_start + (ny - 1) * nx + i
        bot1 = bottom_start + (ny - 1) * nx + i + 1
        add_quad(bot1, bot0, top0, top1)

    for j in range(ny - 1):
        top0 = top_start + j * nx
        top1 = top_start + (j + 1) * nx
        bot0 = bottom_start + j * nx
        bot1 = bottom_start + (j + 1) * nx
        add_quad(bot1, bot0, top0, top1)

        top0 = top_start + j * nx + (nx - 1)
        top1 = top_start + (j + 1) * nx + (nx - 1)
        bot0 = bottom_start + j * nx + (nx - 1)
        bot1 = bottom_start + (j + 1) * nx + (nx - 1)
        add_quad(bot0, bot1, top1, top0)

    return Mesh("terrain_solid_box", vertices, facets)


def building_cell_index(features: list[dict]) -> tuple[list[Polygon], list[float]]:
    polygons: list[Polygon] = []
    heights: list[float] = []
    for feature in features:
        props = feature["properties"]
        mark_z = numeric(props, "MARK_Z")
        top_z = numeric(props, "TAK_Z")
        bygg_h = numeric(props, "BYGG_H")
        if mark_z is None:
            continue
        if top_z is None and bygg_h is not None:
            top_z = mark_z + bygg_h
        if top_z is None or top_z <= mark_z:
            continue
        height = top_z - mark_z
        for polygon in rings_from_geometry(feature["geometry"]):
            if not polygon:
                continue
            outer = clean_ring(polygon[0])
            if len(outer) < 3:
                continue
            poly = Polygon(outer)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if not poly.is_empty and poly.area > 0:
                polygons.append(poly)
                heights.append(height)
    return polygons, heights


def urban_grid_surface_mesh(
    xs: np.ndarray,
    ys: np.ndarray,
    z: np.ndarray,
    building_features: list[dict],
) -> tuple[Mesh, dict]:
    xmin, ymin, _, _ = PILOT_BBOX
    polygons, heights = building_cell_index(building_features)
    tree = STRtree(polygons) if polygons else None
    nx = len(xs)
    ny = len(ys)
    cell_kind: list[list[int | None]] = [[None for _ in range(nx - 1)] for _ in range(ny - 1)]
    cell_height = np.zeros((ny - 1, nx - 1), dtype=float)
    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []

    def terrain_corner(j: int, i: int) -> float:
        return float(z[j, i])

    def add_vertex(x: float, y: float, vz: float) -> int:
        vertices.append((float(x - xmin), float(y - ymin), float(vz)))
        return len(vertices) - 1

    def add_quad_by_points(points: list[tuple[float, float, float]]) -> None:
        base = len(vertices)
        for x, y, vz in points:
            vertices.append((float(x - xmin), float(y - ymin), float(vz)))
        facets.append((base, base + 1, base + 2))
        facets.append((base, base + 2, base + 3))

    # Classify each terrain cell by center point. If multiple buildings overlap,
    # use the tallest; this avoids hidden lower roofs in the grid surface.
    building_cells = 0
    for j in range(ny - 1):
        for i in range(nx - 1):
            cx = float((xs[i] + xs[i + 1]) * 0.5)
            cy = float((ys[j] + ys[j + 1]) * 0.5)
            if tree is not None:
                point = Point(cx, cy)
                candidates = [
                    int(index)
                    for index in tree.query(point)
                    if polygons[int(index)].contains(point)
                ]
                if candidates:
                    best = max(candidates, key=lambda idx: heights[idx])
                    cell_kind[j][i] = best
                    cell_height[j, i] = heights[best]
                    building_cells += 1

    top_edges: dict[tuple[int, int, str], tuple[tuple[float, float, float], tuple[float, float, float]]] = {}

    for j in range(ny - 1):
        for i in range(nx - 1):
            x0 = float(xs[i])
            x1 = float(xs[i + 1])
            y0 = float(ys[j])
            y1 = float(ys[j + 1])
            if cell_kind[j][i] is None:
                p00 = (x0, y0, terrain_corner(j, i))
                p10 = (x1, y0, terrain_corner(j, i + 1))
                p11 = (x1, y1, terrain_corner(j + 1, i + 1))
                p01 = (x0, y1, terrain_corner(j + 1, i))
            else:
                base_z = max(
                    terrain_corner(j, i),
                    terrain_corner(j, i + 1),
                    terrain_corner(j + 1, i + 1),
                    terrain_corner(j + 1, i),
                )
                roof_z = base_z + float(cell_height[j, i])
                p00 = (x0, y0, roof_z)
                p10 = (x1, y0, roof_z)
                p11 = (x1, y1, roof_z)
                p01 = (x0, y1, roof_z)
            add_quad_by_points([p00, p10, p11, p01])
            top_edges[(j, i, "S")] = (p00, p10)
            top_edges[(j, i, "E")] = (p10, p11)
            top_edges[(j, i, "N")] = (p11, p01)
            top_edges[(j, i, "W")] = (p01, p00)

    def add_wall(edge_a, edge_b) -> None:
        a0, a1 = edge_a
        b0, b1 = edge_b
        if max(abs(a0[2] - b1[2]), abs(a1[2] - b0[2])) < 1e-6:
            return
        add_quad_by_points([a0, a1, b0, b1])

    walls = 0
    for j in range(ny - 1):
        for i in range(nx - 1):
            if i + 1 < nx - 1:
                before = len(facets)
                add_wall(top_edges[(j, i, "E")], top_edges[(j, i + 1, "W")])
                walls += (len(facets) - before) // 2
            if j + 1 < ny - 1:
                before = len(facets)
                add_wall(top_edges[(j, i, "N")], top_edges[(j + 1, i, "S")])
                walls += (len(facets) - before) // 2

    return Mesh("urban_grid_surface", vertices, facets), {
        "building_cells": building_cells,
        "vertical_wall_quads": walls,
        "method": "regular DEM cells; building cells are raised to roof height with vertical walls along cell transitions",
    }


def urban_heightfield(
    xs: np.ndarray,
    ys: np.ndarray,
    z: np.ndarray,
    building_features: list[dict],
) -> tuple[np.ndarray, dict]:
    polygons, heights = building_cell_index(building_features)
    tree = STRtree(polygons) if polygons else None
    z_urban = np.array(z, copy=True)
    nx = len(xs)
    ny = len(ys)
    raised_vertices = 0
    if tree is None:
        return z_urban, {"raised_vertices": 0, "method": "terrain only"}

    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            point = Point(float(x), float(y))
            candidates = [
                int(index)
                for index in tree.query(point)
                if polygons[int(index)].buffer(STITCHED_EDGE_M * 0.55).contains(point)
            ]
            if not candidates:
                continue
            best_height = max(heights[index] for index in candidates)
            cell_js = [jj for jj in (j - 1, j) if 0 <= jj < ny - 1]
            cell_is = [ii for ii in (i - 1, i) if 0 <= ii < nx - 1]
            local_ground = max(float(z[jj, ii]) for jj in [j] for ii in [i])
            if cell_js and cell_is:
                local_ground = max(float(z[jj + dj, ii + di]) for jj in cell_js for ii in cell_is for dj in (0, 1) for di in (0, 1))
            z_urban[j, i] = max(z_urban[j, i], local_ground + best_height)
            raised_vertices += 1

    return z_urban, {
        "raised_vertices": raised_vertices,
        "method": "continuous DEM heightfield; vertices near building footprints are raised by building height",
    }


def urban_grid_box_mesh(
    xs: np.ndarray,
    ys: np.ndarray,
    urban_surface: Mesh,
    top_z: float,
) -> Mesh:
    xmin, ymin, _, _ = PILOT_BBOX
    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []

    def add_quad(points: list[tuple[float, float, float]]) -> None:
        base = len(vertices)
        vertices.extend(points)
        facets.append((base, base + 1, base + 2))
        facets.append((base, base + 2, base + 3))

    # Full top cap grid matching the DEM footprint.
    nx = len(xs)
    ny = len(ys)
    for j in range(ny - 1):
        for i in range(nx - 1):
            x0 = float(xs[i] - xmin)
            x1 = float(xs[i + 1] - xmin)
            y0 = float(ys[j] - ymin)
            y1 = float(ys[j + 1] - ymin)
            add_quad([(x0, y0, top_z), (x1, y0, top_z), (x1, y1, top_z), (x0, y1, top_z)])

    # Side walls from urban surface perimeter to the top cap.
    edge_counts: dict[tuple[tuple[float, float, float], tuple[float, float, float]], int] = {}
    edge_oriented: dict[tuple[tuple[float, float, float], tuple[float, float, float]], tuple[tuple[float, float, float], tuple[float, float, float]]] = {}
    for ia, ib, ic in urban_surface.facets:
        tri = [urban_surface.vertices[ia], urban_surface.vertices[ib], urban_surface.vertices[ic]]
        for a, b in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]:
            key = tuple(sorted((a, b)))
            edge_counts[key] = edge_counts.get(key, 0) + 1
            edge_oriented[key] = (a, b)
    for key, count in edge_counts.items():
        if count != 1:
            continue
        a, b = edge_oriented[key]
        ax, ay, _ = a
        bx, by, _ = b
        if (
            abs(ax) < 1e-5 and abs(bx) < 1e-5
            or abs(ay) < 1e-5 and abs(by) < 1e-5
            or abs(ax - (xs[-1] - xmin)) < 1e-5 and abs(bx - (xs[-1] - xmin)) < 1e-5
            or abs(ay - (ys[-1] - ymin)) < 1e-5 and abs(by - (ys[-1] - ymin)) < 1e-5
        ):
            add_quad([a, b, (b[0], b[1], top_z), (a[0], a[1], top_z)])

    return Mesh("urban_grid_box", vertices, facets)


def close_surface_to_top_box(surface: Mesh, top_z: float) -> Mesh:
    edge_counts: dict[tuple[int, int], int] = {}
    edge_oriented: dict[tuple[int, int], tuple[int, int]] = {}
    for ia, ib, ic in surface.facets:
        for a, b in ((ia, ib), (ib, ic), (ic, ia)):
            key = tuple(sorted((a, b)))
            edge_counts[key] = edge_counts.get(key, 0) + 1
            edge_oriented[key] = (a, b)

    boundary_edges = [edge_oriented[key] for key, count in edge_counts.items() if count == 1]
    if not boundary_edges:
        return Mesh("top_box_closure", [], [])

    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []

    def add_vertex(point: tuple[float, float, float]) -> int:
        vertices.append(point)
        return len(vertices) - 1

    top_vertex: dict[tuple[float, float], int] = {}

    def get_top(x: float, y: float) -> int:
        key = (round(x, 6), round(y, 6))
        if key not in top_vertex:
            top_vertex[key] = add_vertex((x, y, top_z))
        return top_vertex[key]

    # Side walls along all boundary loops.
    for ia, ib in boundary_edges:
        a = surface.vertices[ia]
        b = surface.vertices[ib]
        va = add_vertex(a)
        vb = add_vertex(b)
        tb = get_top(b[0], b[1])
        ta = get_top(a[0], a[1])
        facets.append((va, vb, tb))
        facets.append((va, tb, ta))

    # Top cap: triangulate the x/y footprint of the boundary loops.
    xy_points = sorted(top_vertex.keys())
    if len(xy_points) >= 3:
        polygon = Polygon(xy_points).convex_hull
        top_tris = shapely.constrained_delaunay_triangles(polygon)
        for tri in top_tris.geoms:
            coords = list(tri.exterior.coords)[:-1]
            if len(coords) != 3:
                continue
            a = get_top(float(coords[0][0]), float(coords[0][1]))
            b = get_top(float(coords[1][0]), float(coords[1][1]))
            c = get_top(float(coords[2][0]), float(coords[2][1]))
            facets.append((a, c, b))

    return Mesh("top_box_closure", vertices, facets)


def numeric(properties: dict, key: str) -> float | None:
    value = properties.get(key)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def draped_base_z(ring: list[tuple[float, float]], sampler: TerrainSampler) -> float:
    samples = [sampler.sample(x, y) for x, y in ring]
    cx = sum(x for x, _ in ring) / len(ring)
    cy = sum(y for _, y in ring) / len(ring)
    samples.append(sampler.sample(cx, cy))
    return max(samples)


def add_building_part(
    ring: list[tuple[float, float]],
    bottom_z: float,
    top_z: float,
    vertices: list[tuple[float, float, float]],
    facets: list[tuple[int, int, int]],
) -> bool:
    xmin, ymin, _, _ = PILOT_BBOX
    ring = clean_ring(ring)
    if len(ring) < 3:
        return True
    if polygon_area(ring) < 0:
        ring.reverse()

    base = len(vertices)
    for x, y in ring:
        vertices.append((x - xmin, y - ymin, bottom_z))
    top_base = len(vertices)
    for x, y in ring:
        vertices.append((x - xmin, y - ymin, top_z))

    triangles, clean = triangulate_ring(ring)
    for a, b, c in triangles:
        facets.append((top_base + a, top_base + b, top_base + c))
        facets.append((base + c, base + b, base + a))

    n = len(ring)
    for i in range(n):
        j = (i + 1) % n
        facets.append((base + i, base + j, top_base + j))
        facets.append((base + i, top_base + j, top_base + i))
    return clean


def add_building_surface_part(
    ring: list[tuple[float, float]],
    height: float,
    terrain_sampler: TerrainSampler,
    vertices: list[tuple[float, float, float]],
    facets: list[tuple[int, int, int]],
) -> bool:
    xmin, ymin, _, _ = PILOT_BBOX
    ring = densify_ring(ring)
    if len(ring) < 3:
        return True
    if polygon_area(ring) < 0:
        ring.reverse()

    base = len(vertices)
    bottom_z: list[float] = []
    for x, y in ring:
        z = terrain_sampler.sample(x, y)
        bottom_z.append(z)
        vertices.append((x - xmin, y - ymin, z))
    top_base = len(vertices)
    for (x, y), z in zip(ring, bottom_z):
        vertices.append((x - xmin, y - ymin, z + height))

    triangles, clean = triangulate_ring(ring)
    for a, b, c in triangles:
        facets.append((top_base + a, top_base + b, top_base + c))

    n = len(ring)
    for i in range(n):
        j = (i + 1) % n
        facets.append((base + i, base + j, top_base + j))
        facets.append((base + i, top_base + j, top_base + i))
    return clean


def buildings_mesh(features: list[dict], terrain_sampler: TerrainSampler | None = None) -> tuple[Mesh, dict]:
    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []
    skipped = 0
    holes = 0
    fallback_triangulations = 0
    draped = 0
    base_deltas: list[float] = []
    suspicious: list[dict] = []

    for idx, feature in enumerate(features, start=1):
        props = feature["properties"]
        mark_z = numeric(props, "MARK_Z")
        top_z = numeric(props, "TAK_Z")
        bygg_h = numeric(props, "BYGG_H")
        if mark_z is None:
            skipped += 1
            continue
        if top_z is None and bygg_h is not None:
            top_z = mark_z + bygg_h
        if top_z is None or top_z <= mark_z:
            skipped += 1
            continue
        height = top_z - mark_z
        if bygg_h is not None and bygg_h >= SUSPICIOUS_HEIGHT_M:
            suspicious.append(
                {
                    "feature_index": idx,
                    "UUID": props.get("UUID"),
                    "GRUPP": props.get("GRUPP"),
                    "BYGG_H": bygg_h,
                    "MARK_Z": mark_z,
                    "TAK_Z": top_z,
                }
            )

        for polygon in rings_from_geometry(feature["geometry"]):
            if not polygon:
                continue
            holes += max(0, len(polygon) - 1)
            ring = clean_ring(polygon[0])
            bottom_z = mark_z
            if terrain_sampler is not None and len(ring) >= 3:
                bottom_z = draped_base_z(ring, terrain_sampler)
                draped += 1
                base_deltas.append(bottom_z - mark_z)
            clean = add_building_part(ring, bottom_z, bottom_z + height, vertices, facets)
            if not clean:
                fallback_triangulations += 1

    return Mesh("buildings", vertices, facets), {
        "skipped_buildings": skipped,
        "holes_ignored": holes,
        "fallback_triangulations": fallback_triangulations,
        "draped_building_parts": draped,
        "base_delta_min": min(base_deltas) if base_deltas else None,
        "base_delta_mean": sum(base_deltas) / len(base_deltas) if base_deltas else None,
        "base_delta_max": max(base_deltas) if base_deltas else None,
        "suspicious_buildings": suspicious,
    }


def stitched_buildings_mesh(features: list[dict], terrain_sampler: TerrainSampler) -> tuple[Mesh, dict]:
    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []
    skipped = 0
    holes = 0
    fallback_triangulations = 0
    parts = 0

    for feature in features:
        props = feature["properties"]
        mark_z = numeric(props, "MARK_Z")
        top_z = numeric(props, "TAK_Z")
        bygg_h = numeric(props, "BYGG_H")
        if mark_z is None:
            skipped += 1
            continue
        if top_z is None and bygg_h is not None:
            top_z = mark_z + bygg_h
        if top_z is None or top_z <= mark_z:
            skipped += 1
            continue
        height = top_z - mark_z
        for polygon in rings_from_geometry(feature["geometry"]):
            if not polygon:
                continue
            holes += max(0, len(polygon) - 1)
            clean = add_building_surface_part(polygon[0], height, terrain_sampler, vertices, facets)
            parts += 1
            if not clean:
                fallback_triangulations += 1

    return Mesh("stitched_buildings", vertices, facets), {
        "skipped_buildings": skipped,
        "holes_ignored": holes,
        "fallback_triangulations": fallback_triangulations,
        "stitched_building_parts": parts,
    }


def stitched_building_regions_mesh(
    regions: list[Polygon],
    heights: list[float],
    terrain_sampler: TerrainSampler,
) -> tuple[Mesh, dict]:
    vertices: list[tuple[float, float, float]] = []
    facets: list[tuple[int, int, int]] = []
    fallback_triangulations = 0
    parts = 0
    for region, height in zip(regions, heights):
        clean = add_building_surface_part(
            list(region.exterior.coords)[:-1],
            height,
            terrain_sampler,
            vertices,
            facets,
        )
        parts += 1
        if not clean:
            fallback_triangulations += 1
        for interior in region.interiors:
            clean = add_building_surface_part(
                list(interior.coords)[:-1],
                height,
                terrain_sampler,
                vertices,
                facets,
            )
            parts += 1
            if not clean:
                fallback_triangulations += 1
    return Mesh("stitched_building_regions", vertices, facets), {
        "region_parts": parts,
        "fallback_triangulations": fallback_triangulations,
    }


def facet_normal(a: tuple[float, float, float], b: tuple[float, float, float], c: tuple[float, float, float]) -> tuple[float, float, float]:
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    length = math.sqrt(nx * nx + ny * ny + nz * nz)
    if length == 0:
        return (0.0, 0.0, 0.0)
    return (nx / length, ny / length, nz / length)


def write_binary_stl(path: Path, meshes: list[Mesh]) -> int:
    facet_count = sum(len(mesh.facets) for mesh in meshes)
    with path.open("wb") as f:
        header = f"Urban pilot STL: {path.name}".encode("ascii", errors="ignore")[:80]
        f.write(header.ljust(80, b" "))
        f.write(struct.pack("<I", facet_count))
        for mesh in meshes:
            for ia, ib, ic in mesh.facets:
                a = mesh.vertices[ia]
                b = mesh.vertices[ib]
                c = mesh.vertices[ic]
                normal = facet_normal(a, b, c)
                f.write(struct.pack("<12fH", *normal, *a, *b, *c, 0))
    return facet_count


def render_preview(
    xs: np.ndarray,
    ys: np.ndarray,
    z: np.ndarray,
    building_features: list[dict],
    out_path: Path,
    terrain_sampler: TerrainSampler | None = None,
) -> None:
    xmin, ymin, _, _ = PILOT_BBOX
    xx, yy = np.meshgrid(xs - xmin, ys - ymin)

    fig = plt.figure(figsize=(12, 9), dpi=180)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(xx, yy, z, cmap="terrain", linewidth=0, antialiased=True, alpha=0.82, rstride=3, cstride=3)

    for feature in building_features:
        props = feature["properties"]
        mark_z = numeric(props, "MARK_Z")
        top_z = numeric(props, "TAK_Z")
        if mark_z is None or top_z is None:
            continue
        height = top_z - mark_z
        for polygon in rings_from_geometry(feature["geometry"]):
            if not polygon:
                continue
            ring = clean_ring(polygon[0])
            if len(ring) < 3:
                continue
            draw_top_z = top_z
            if terrain_sampler is not None:
                draw_top_z = draped_base_z(ring, terrain_sampler) + height
            x = [p[0] - xmin for p in ring] + [ring[0][0] - xmin]
            y = [p[1] - ymin for p in ring] + [ring[0][1] - ymin]
            ax.plot(x, y, [draw_top_z] * len(x), color="#0f766e", linewidth=0.45, alpha=0.9)
            if height >= 20:
                cx = sum(x[:-1]) / (len(x) - 1)
                cy = sum(y[:-1]) / (len(y) - 1)
                ax.plot([cx, cx], [cy, cy], [draw_top_z - height, draw_top_z], color="#7f1d1d", linewidth=0.7, alpha=0.8)

    ax.set_title(f"{OUTPUT_LABEL} STL preview: terrain and LOD1 building roofs")
    ax.set_xlabel("local x (m)")
    ax.set_ylabel("local y (m)")
    ax.set_zlabel("elevation (m)")
    ax.view_init(elev=38, azim=-132)
    ax.set_box_aspect((1, 1, 0.16))
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def write_metadata(path: Path, metadata: dict) -> None:
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    require_tools(["ogr2ogr", "gdal_grid", "gdal_translate", "gmsh"])
    out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    buildings_geojson = out_dir / "buildings_epsg3006.geojson"
    terrain_geojson = out_dir / "terrain_contours_epsg3006.geojson"
    export_layer(INPUT_GPKG, BUILDING_LAYER, buildings_geojson)
    export_layer(INPUT_GPKG, TERRAIN_LAYER, terrain_geojson)

    building_features = load_geojson(buildings_geojson)
    terrain_features = load_geojson(terrain_geojson)
    points = densified_contour_points(terrain_features)
    vrt_path = write_contour_vrt(points, out_dir)
    _, xyz_path = grid_terrain(vrt_path, out_dir)
    xs, ys, z = load_dem_xyz(xyz_path)

    t_mesh = terrain_mesh(xs, ys, z)
    terrain_sampler = TerrainSampler(xs, ys, z)
    b_mesh, building_meta = buildings_mesh(building_features, terrain_sampler)
    union_regions, union_heights, union_meta = building_union_regions(building_features)
    cgal_regions, cgal_heights, cgal_filter_meta = filter_interior_building_regions(union_regions, union_heights)

    terrain_box_mesh_obj = terrain_solid_box_mesh(
        xs, ys, z, float(np.nanmin(z)) - TERRAIN_BOTTOM_MARGIN_M
    )
    terrain_box_facets = write_binary_stl(out_dir / "terrain_box.stl", [terrain_box_mesh_obj])
    combined_stl = out_dir / f"{OUTPUT_SLUG}_combined.stl"
    boxed_stl = out_dir / f"{OUTPUT_SLUG}_boxed.stl"
    preview_png = out_dir / f"{OUTPUT_SLUG}_preview.png"
    cgal_stats = run_cgal_stl(out_dir, xs, ys, z, cgal_regions, cgal_heights, BOX_TOP_Z)
    render_preview(xs, ys, z, building_features, preview_png, terrain_sampler)

    heights = [numeric(f["properties"], "BYGG_H") for f in building_features]
    heights = [h for h in heights if h is not None]
    contour_z = [numeric(f["properties"], "HOJD") for f in terrain_features]
    contour_z = [h for h in contour_z if h is not None]
    metadata = {
        "label": OUTPUT_LABEL,
        "source_gpkg": str(INPUT_GPKG),
        "input_layers": {
            "buildings": BUILDING_LAYER,
            "terrain_contours": TERRAIN_LAYER,
        },
        "bbox_epsg3006": PILOT_BBOX,
        "bbox_size_m": {
            "x": PILOT_BBOX[2] - PILOT_BBOX[0],
            "y": PILOT_BBOX[3] - PILOT_BBOX[1],
        },
        "local_origin_epsg3006": {"x": PILOT_BBOX[0], "y": PILOT_BBOX[1]},
        "grid_size": [len(xs), len(ys)],
        "terrain_grid_spacing_m": {
            "x": float(abs(xs[1] - xs[0])) if len(xs) > 1 else None,
            "y": float(abs(ys[1] - ys[0])) if len(ys) > 1 else None,
        },
        "terrain_contours": {
            "features": len(terrain_features),
            "densified_points": len(points),
            "min_z": min(contour_z),
            "max_z": max(contour_z),
        },
        "terrain_dem_z": {
            "min": float(np.nanmin(z)),
            "max": float(np.nanmax(z)),
            "nan_count": int(np.isnan(z).sum()),
        },
        "buildings": {
            "features": len(building_features),
            "height_min": min(heights),
            "height_mean": sum(heights) / len(heights),
            "height_max": max(heights),
            **building_meta,
        },
        "stl_facets": {
            "terrain": cgal_stats["terrain_facets"],
            "terrain_box": terrain_box_facets,
            "buildings": cgal_stats["buildings_facets"],
            "combined": cgal_stats["combined_facets"],
            "boxed": cgal_stats["boxed_facets"],
        },
        "boxed": {
            "top_z": BOX_TOP_Z,
            "description": "CGAL constrained urban surface plus vertical side walls and a flat top cap.",
            "uses_carved_terrain": False,
        },
        "terrain_box": {
            "bottom_z": float(np.nanmin(z)) - TERRAIN_BOTTOM_MARGIN_M,
            "description": "Closed terrain solid with interpolated terrain top, flat bottom cap, and vertical side walls.",
        },
        "cgal_constrained_combined": {
            "enabled": True,
            "union_regions": union_meta,
            "region_filter": cgal_filter_meta,
            "cdt_vertices": cgal_stats["cdt_vertices"],
            "cdt_faces": cgal_stats["cdt_faces"],
            "building_regions": cgal_stats["building_regions"],
            "method": "CGAL Constrained_Delaunay_triangulation_2 with building footprints as hard constraints and holes in the terrain domain.",
        },
        "outputs": {
            "terrain_stl": str(out_dir / "terrain.stl"),
            "terrain_box_stl": str(out_dir / "terrain_box.stl"),
            "buildings_stl": str(out_dir / "buildings.stl"),
            "combined_stl": str(combined_stl),
            "boxed_stl": str(boxed_stl),
            "preview_png": str(preview_png),
        },
    }
    write_metadata(out_dir / "metadata.json", metadata)

    print("\nWrote:")
    for key, value in metadata["outputs"].items():
        print(f"  {key}: {value}")
    print(f"  metadata: {out_dir / 'metadata.json'}")
    print(f"\nCombined facets: {cgal_stats['combined_facets']:,}")
    if building_meta["suspicious_buildings"]:
        print("Suspicious buildings retained:")
        for item in building_meta["suspicious_buildings"]:
            print(f"  {item}")


if __name__ == "__main__":
    main()
