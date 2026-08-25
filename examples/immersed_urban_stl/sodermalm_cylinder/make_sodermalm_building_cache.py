#!/usr/bin/env python3
"""Write a Brinkman indicator cache from Sodermalm building footprints."""

import argparse
import json
import math
import struct
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


HEADER_SIZE = 132
MARKER_SIZE = 4
SHORE_BLEND_M = 180.0
BUILDING_MIN_TERRAIN_ABOVE_WATER_M = 0.75
BUILDING_BASE_CLEARANCE_M = 0.5
BUILDING_MIN_FOOTPRINT_AREA_M2 = 0.0
BUILDING_MIN_EDGE_M = 0.0
BUILDING_SIMPLIFY_M = 0.0
BUILDING_MIN_HEIGHT_M = 2.0
BUILDING_MAX_HEIGHT_M = 65.0


def load_geojson(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def line_points(geometry: Dict[str, object]) -> List[Tuple[float, float]]:
    if not geometry:
        return []
    if geometry["type"] == "LineString":
        return [(float(row[0]), float(row[1])) for row in geometry["coordinates"]]
    if geometry["type"] == "MultiLineString":
        return [
            (float(row[0]), float(row[1]))
            for part in geometry["coordinates"]
            for row in part
        ]
    return []


def polygon_rings(geometry: Dict[str, object]) -> List[List[Tuple[float, float]]]:
    if not geometry:
        return []
    if geometry["type"] == "Polygon":
        return [[(float(row[0]), float(row[1])) for row in geometry["coordinates"][0]]]
    if geometry["type"] == "MultiPolygon":
        return [
            [(float(row[0]), float(row[1])) for row in poly[0]]
            for poly in geometry["coordinates"]
        ]
    return []


def polygon_boundary_rings(geometry: Dict[str, object]) -> List[List[Tuple[float, float]]]:
    if not geometry:
        return []
    if geometry["type"] == "Polygon":
        return [
            [(float(row[0]), float(row[1])) for row in ring]
            for ring in geometry["coordinates"]
        ]
    if geometry["type"] == "MultiPolygon":
        return [
            [(float(row[0]), float(row[1])) for row in ring]
            for poly in geometry["coordinates"]
            for ring in poly
        ]
    return []


def point_in_polygon(x: float, y: float, ring: List[Tuple[float, float]]) -> bool:
    inside = False
    for i, p0 in enumerate(ring):
        x0, y0 = p0
        x1, y1 = ring[(i + 1) % len(ring)]
        if (y0 > y) != (y1 > y):
            x_cross = (x1 - x0) * (y - y0) / (y1 - y0 + 1.0e-300) + x0
            if x < x_cross:
                inside = not inside
    return inside


def point_in_land_mask(x: float, y: float, boundary_rings: List[List[Tuple[float, float]]]) -> bool:
    inside = False
    for ring in boundary_rings:
        if point_in_polygon(x, y, ring):
            inside = not inside
    return inside


def point_segment_distance(px: float, py: float, ax: float, ay: float, bx: float, by: float) -> float:
    vx = bx - ax
    vy = by - ay
    wx = px - ax
    wy = py - ay
    denom = vx * vx + vy * vy
    t = 0.0 if denom == 0.0 else max(0.0, min(1.0, (wx * vx + wy * vy) / denom))
    qx = ax + t * vx
    qy = ay + t * vy
    return math.hypot(px - qx, py - qy)


def distance_to_ring(x: float, y: float, ring: List[Tuple[float, float]]) -> float:
    return min(
        point_segment_distance(x, y, ring[i][0], ring[i][1], ring[(i + 1) % len(ring)][0], ring[(i + 1) % len(ring)][1])
        for i in range(len(ring) - 1)
    )


def smoothstep(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


def terrain_samples(contours: Path, center: Tuple[float, float]) -> Tuple[np.ndarray, float]:
    rows = []
    for feature in load_geojson(contours)["features"]:
        z = feature["properties"].get("HOJD")
        if z is None:
            continue
        for x, y in line_points(feature["geometry"]):
            rows.append((x - center[0], y - center[1], float(z)))
    if len(rows) < 10:
        raise RuntimeError("Too few terrain contour samples for Sodermalm cylinder.")
    data = np.array(rows, dtype=float)
    return data, float(np.min(data[:, 2]))


def idw_height(x: float, y: float, samples: np.ndarray) -> float:
    xy = samples[:, :2]
    zz = samples[:, 2]
    dist2 = (xy[:, 0] - x) ** 2 + (xy[:, 1] - y) ** 2
    idx = np.argpartition(dist2, min(20, len(dist2) - 1))[:20]
    d = np.sqrt(dist2[idx])
    if float(d.min()) < 1.0e-9:
        return float(zz[idx[int(d.argmin())]])
    w = 1.0 / np.maximum(d, 2.0) ** 2
    return float(np.sum(w * zz[idx]) / np.sum(w))


def terrain_height(
    x: float,
    y: float,
    samples: np.ndarray,
    water_level: float,
    shoreline: List[List[Tuple[float, float]]],
    center: Tuple[float, float],
    radius: float,
) -> float:
    gx = x + center[0]
    gy = y + center[1]
    if not point_in_land_mask(gx, gy, shoreline):
        return water_level
    d = min(distance_to_ring(gx, gy, ring) for ring in shoreline)
    weight = smoothstep(d / SHORE_BLEND_M)
    raw = max(water_level, idw_height(x, y, samples))
    edge_weight = smoothstep(max(0.0, radius - math.hypot(x, y)) / SHORE_BLEND_M)
    return water_level + weight * edge_weight * (raw - water_level)


def ring_centroid(ring: List[Tuple[float, float]]) -> Tuple[float, float]:
    if not ring:
        return (0.0, 0.0)
    pts = ring[:-1] if ring[0] == ring[-1] else ring
    area = 0.0
    cx = 0.0
    cy = 0.0
    for p0, p1 in zip(pts, pts[1:] + pts[:1]):
        x0, y0 = p0
        x1, y1 = p1
        cross = x0 * y1 - x1 * y0
        area += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    if abs(area) < 1.0e-9:
        return (sum(x for x, _ in pts) / len(pts), sum(y for _, y in pts) / len(pts))
    return (cx / (3.0 * area), cy / (3.0 * area))


def clean_building_ring(ring: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    cleaned = []
    for point in ring:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def building_ring_area(ring: List[Tuple[float, float]]) -> float:
    return 0.5 * sum(
        p0[0] * p1[1] - p1[0] * p0[1]
        for p0, p1 in zip(ring, ring[1:] + ring[:1])
    )


def building_ring_min_edge(ring: List[Tuple[float, float]]) -> float:
    if len(ring) < 2:
        return 0.0
    return min(
        math.hypot(p1[0] - p0[0], p1[1] - p0[1])
        for p0, p1 in zip(ring, ring[1:] + ring[:1])
    )


def perpendicular_distance(point: Tuple[float, float], start: Tuple[float, float], end: Tuple[float, float]) -> float:
    return point_segment_distance(point[0], point[1], start[0], start[1], end[0], end[1])


def simplify_ring(points: List[Tuple[float, float]], tolerance: float) -> List[Tuple[float, float]]:
    pts = points[:-1] if points and points[0] == points[-1] else list(points)
    if len(pts) <= 4 or tolerance <= 0.0:
        return pts

    def rdp(rows: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
        if len(rows) <= 2:
            return rows
        start = rows[0]
        end = rows[-1]
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


def building_base_and_top(properties: Dict[str, object], terrain_base: float, water_level: float) -> Tuple[float, float]:
    mark = properties.get("MARK_Z")
    roof = properties.get("TAK_Z")
    by_height = properties.get("BYGG_H")
    base_z = max(water_level, terrain_base + BUILDING_BASE_CLEARANCE_M)

    try:
        height = float(roof) - float(mark)
        if BUILDING_MIN_HEIGHT_M <= height <= BUILDING_MAX_HEIGHT_M:
            return base_z, base_z + height
    except (TypeError, ValueError):
        pass

    try:
        height = float(by_height)
        if BUILDING_MIN_HEIGHT_M <= height <= BUILDING_MAX_HEIGHT_M:
            return base_z, base_z + height
    except (TypeError, ValueError):
        pass
    return base_z, base_z + 12.0


def building_sits_on_water_level_terrain(terrain_base: float, water_level: float) -> bool:
    return terrain_base <= water_level + BUILDING_MIN_TERRAIN_ABOVE_WATER_M


def parse_fld_header(path: Path) -> Dict[str, object]:
    with path.open("rb") as handle:
        raw = handle.read(HEADER_SIZE)
        marker = handle.read(MARKER_SIZE)
    header = raw.decode("ascii", errors="replace")
    parts = header.split()
    if len(parts) < 11 or parts[0] != "#std":
        raise RuntimeError(f"Unexpected Neko field header in {path}: {header!r}")
    little = round(struct.unpack("<f", marker)[0], 5)
    big = round(struct.unpack(">f", marker)[0], 5)
    if little == 6.54321:
        endian = "<"
    elif big == 6.54321:
        endian = ">"
    else:
        raise RuntimeError(f"Could not determine Neko field endianness for {path}")
    return {
        "wdsz": int(parts[1]),
        "lx": int(parts[2]),
        "ly": int(parts[3]),
        "lz": int(parts[4]),
        "nelv": int(parts[5]),
        "nelgt": int(parts[6]),
        "time": float(parts[7].replace("D", "E")),
        "step": int(parts[8]),
        "fid0": int(parts[9]),
        "nfileo": int(parts[10]),
        "rdcode": header[83:93],
        "endian": endian,
    }


def field_offsets(header: Dict[str, object]) -> Dict[str, int]:
    npts = header["lx"] * header["ly"] * header["lz"]
    nelv = header["nelv"]
    scalar_bytes = npts * header["wdsz"]
    offset = HEADER_SIZE + MARKER_SIZE + 4 * nelv
    offsets = {}
    rdcode = header["rdcode"].strip()
    i = 0
    while i < len(rdcode):
        code = rdcode[i]
        if code == "X":
            offsets["X"] = offset
            offset += nelv * 3 * scalar_bytes
        elif code == "U":
            offsets["U"] = offset
            offset += nelv * 3 * scalar_bytes
        elif code in ("P", "T"):
            offsets[code] = offset
            offset += nelv * scalar_bytes
        elif code == "S":
            n_scalars = int(rdcode[i + 1 : i + 3])
            offsets["S"] = offset
            offset += nelv * n_scalars * scalar_bytes
            i += 2
        i += 1
    return offsets


def read_template_coordinates(path: Path) -> Tuple[Dict[str, object], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    header = parse_fld_header(path)
    offsets = field_offsets(header)
    if "X" not in offsets:
        raise RuntimeError(f"Template field {path} does not contain mesh coordinates")

    npts = header["lx"] * header["ly"] * header["lz"]
    count = header["nelv"] * npts
    dtype = np.dtype(f"{header['endian']}f{header['wdsz']}")
    with path.open("rb") as handle:
        handle.seek(HEADER_SIZE + MARKER_SIZE)
        element_ids = np.fromfile(handle, dtype=np.dtype(f"{header['endian']}i4"), count=header["nelv"])
        handle.seek(offsets["X"])
        coords = np.fromfile(handle, dtype=dtype, count=3 * count).astype(np.float64, copy=False)
    if element_ids.size != header["nelv"] or coords.size != 3 * count:
        raise RuntimeError(f"Template field {path} is truncated")
    coords = coords.reshape(3, count)
    return header, element_ids, coords[0], coords[1], coords[2]


def point_segment_distance_xy(
    x: np.ndarray,
    y: np.ndarray,
    ax: float,
    ay: float,
    bx: float,
    by: float,
) -> np.ndarray:
    vx = bx - ax
    vy = by - ay
    denom = vx * vx + vy * vy
    if denom <= 0.0:
        return np.hypot(x - ax, y - ay)
    t = np.clip(((x - ax) * vx + (y - ay) * vy) / denom, 0.0, 1.0)
    px = ax + t * vx
    py = ay + t * vy
    return np.hypot(x - px, y - py)


def inside_polygon_xy(x: np.ndarray, y: np.ndarray, ring: List[Tuple[float, float]]) -> np.ndarray:
    inside = np.zeros(x.shape, dtype=bool)
    xj, yj = ring[-1]
    for xi, yi in ring:
        crosses = ((yi > y) != (yj > y)) & (
            x < (xj - xi) * (y - yi) / (yj - yi + 1.0e-300) + xi
        )
        inside ^= crosses
        xj, yj = xi, yi
    return inside


def distance_to_polygon_boundary_xy(
    x: np.ndarray,
    y: np.ndarray,
    ring: List[Tuple[float, float]],
) -> np.ndarray:
    distance = np.full(x.shape, np.inf, dtype=np.float64)
    ax, ay = ring[-1]
    for bx, by in ring:
        distance = np.minimum(distance, point_segment_distance_xy(x, y, ax, ay, bx, by))
        ax, ay = bx, by
    return distance


def smooth_step(distance: np.ndarray, width: float) -> np.ndarray:
    if width <= 0.0:
        return (distance <= 0.0).astype(np.float64)
    t = np.clip((distance - width) / -width, 0.0, 1.0)
    return t**3 * (t * (6.0 * t - 15.0) + 10.0)


def building_parts(generated: Path) -> Tuple[List[Dict[str, object]], float]:
    metadata_candidates = sorted(generated.glob("*_metadata.json"))
    if not metadata_candidates:
        raise RuntimeError(f"No metadata JSON found in {generated}")
    metadata = json.loads(metadata_candidates[0].read_text(encoding="utf-8"))
    center = tuple(float(v) for v in metadata["center_epsg3006"])
    radius = float(metadata["radius_m"])

    contours = generated / "sodermalm_cylinder_contours.geojson"
    buildings = generated / "sodermalm_buildings.geojson"
    land = generated / "sodermalm_osm_island_epsg3006.geojson"
    if not contours.exists() or not buildings.exists() or not land.exists():
        raise RuntimeError(f"Missing prepared geometry files in {generated}")

    samples, water_level = terrain_samples(contours, center)
    shoreline = polygon_boundary_rings(load_geojson(land)["features"][0]["geometry"])

    parts = []
    skipped = 0
    skipped_small = 0
    skipped_skinny = 0
    for feature in load_geojson(buildings)["features"]:
        props = feature.get("properties", {})
        for ring_global in polygon_rings(feature["geometry"]):
            pts_global = clean_building_ring(ring_global)
            if len(pts_global) < 3:
                skipped += 1
                continue
            pts_global = clean_building_ring(simplify_ring(pts_global, BUILDING_SIMPLIFY_M))
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
            ring = [(x - center[0], y - center[1]) for x, y in pts_global]
            if building_ring_area(ring) < 0.0:
                ring.reverse()
            xs = [p[0] for p in ring]
            ys = [p[1] for p in ring]
            parts.append(
                {
                    "ring": ring,
                    "base_z": base_z,
                    "top_z": top_z,
                    "bounds": (min(xs), min(ys), max(xs), max(ys)),
                }
            )
    print(
        "Building filter: "
        f"kept={len(parts)} skipped={skipped} "
        f"small={skipped_small} skinny={skipped_skinny} "
        f"min_area={BUILDING_MIN_FOOTPRINT_AREA_M2} min_edge={BUILDING_MIN_EDGE_M} simplify={BUILDING_SIMPLIFY_M}"
    )
    return parts, water_level


def compute_indicator(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    parts: List[Dict[str, object]],
    smooth_width: float,
) -> Tuple[np.ndarray, int]:
    indicator = np.zeros(x.shape, dtype=np.float64)
    touched = 0
    for part in parts:
        xmin, ymin, xmax, ymax = part["bounds"]
        in_box = (
            (x >= xmin - smooth_width)
            & (x <= xmax + smooth_width)
            & (y >= ymin - smooth_width)
            & (y <= ymax + smooth_width)
            & (z >= part["base_z"] - smooth_width)
            & (z <= part["top_z"] + smooth_width)
        )
        if not np.any(in_box):
            continue
        ids = np.flatnonzero(in_box)
        xx = x[ids]
        yy = y[ids]
        zz = z[ids]
        inside_xy = inside_polygon_xy(xx, yy, part["ring"])
        boundary_distance = distance_to_polygon_boundary_xy(xx, yy, part["ring"])
        horizontal_signed = np.where(inside_xy, -boundary_distance, boundary_distance)
        vertical_distance = np.maximum(part["base_z"] - zz, zz - part["top_z"])
        signed_distance = np.maximum(horizontal_signed, vertical_distance)
        values = smooth_step(signed_distance, smooth_width)
        changed = values > indicator[ids]
        if np.any(changed):
            indicator[ids[changed]] = values[changed]
            touched += int(np.count_nonzero(changed))
    return indicator, touched


def write_cache(prefix: Path, header: Dict[str, object], element_ids: np.ndarray, x: np.ndarray, y: np.ndarray, z: np.ndarray, s01: np.ndarray) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    data_path = prefix.parent / f"{prefix.name}0.f00000"
    index_path = prefix.parent / f"{prefix.name}0.nek5000"
    lx, ly, lz = header["lx"], header["ly"], header["lz"]
    nelv = header["nelv"]
    nelgt = header["nelgt"]
    header_text = (
        f"#std 8 {lx:2d} {ly:2d} {lz:2d} {nelv:10d} {nelgt:10d} "
        f"{0.0:20.13E} {0:9d} {1:6d} {1:6d} {'XS01':<10s}"
    )
    header_bytes = header_text.encode("ascii")
    if len(header_bytes) > HEADER_SIZE:
        raise RuntimeError("Generated Neko field header is too long")
    with data_path.open("wb") as handle:
        handle.write(header_bytes.ljust(HEADER_SIZE, b" "))
        handle.write(struct.pack("<f", 6.54321))
        element_ids.astype("<i4", copy=False).tofile(handle)
        for values in (x, y, z, s01):
            values.astype("<f8", copy=False).tofile(handle)
    index_path.write_text(
        f"filetemplate:         {prefix.name}%01d.f%05d\n"
        "firsttimestep:     0\n"
        "numtimesteps:     1\n",
        encoding="ascii",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated", type=Path, required=True)
    parser.add_argument("--template-field", type=Path, required=True)
    parser.add_argument("--cache-prefix", type=Path, required=True)
    parser.add_argument("--smooth-width", type=float, default=50.0)
    parser.add_argument("--min-footprint-area", type=float, default=0.0)
    parser.add_argument("--min-edge", type=float, default=0.0)
    parser.add_argument("--simplify", type=float, default=0.0)
    args = parser.parse_args()

    global BUILDING_MIN_FOOTPRINT_AREA_M2, BUILDING_MIN_EDGE_M, BUILDING_SIMPLIFY_M
    BUILDING_MIN_FOOTPRINT_AREA_M2 = args.min_footprint_area
    BUILDING_MIN_EDGE_M = args.min_edge
    BUILDING_SIMPLIFY_M = args.simplify

    header, element_ids, x, y, z = read_template_coordinates(args.template_field)
    parts, _ = building_parts(args.generated)
    indicator, touched = compute_indicator(x, y, z, parts, args.smooth_width)
    write_cache(args.cache_prefix, header, element_ids, x, y, z, indicator)

    solid = int(np.count_nonzero(indicator >= 0.5))
    print(f"Building parts: {len(parts)}")
    print(f"Candidate updates: {touched}")
    print(f"Solid fraction >= 0.5: {solid / indicator.size:.6f}")
    print(f"Wrote {args.cache_prefix}0.f00000")


if __name__ == "__main__":
    main()
