#!/usr/bin/env python3
"""Check building-base contact against the generated terrain surface.

The building STL uses flat LOD1 solids. This audit reports where those flat
bases float above, or cut into, the terrain-following mesh surface.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

import build_sodermalm_cylinder as base


HERE = Path(__file__).resolve().parent
DEFAULT_GENERATED = HERE / "generated_p7_xy35_inlet220_nz10"
PAD = 34
CANVAS_N = 1600


def metadata_path(generated: Path) -> Path:
    matches = sorted(p for p in generated.glob("*metadata.json") if not p.name.startswith("._"))
    if not matches:
        raise FileNotFoundError(f"No metadata JSON found in {generated}")
    return matches[0]


def load_metadata(generated: Path) -> dict:
    return json.loads(metadata_path(generated).read_text(encoding="utf-8"))


def set_building_constants(meta: dict) -> None:
    stl = meta.get("building_stl", {})
    base.BUILDING_BASE_CLEARANCE_M = float(stl.get("base_clearance_m", base.BUILDING_BASE_CLEARANCE_M))
    base.BUILDING_BASE_EMBED_M = float(stl.get("base_embed_m", base.BUILDING_BASE_EMBED_M))
    base.BUILDING_CONTACT_SAMPLE_SPACING_M = float(
        stl.get("contact_sample_spacing_m", base.BUILDING_CONTACT_SAMPLE_SPACING_M)
    )
    base.BUILDING_MIN_FOOTPRINT_AREA_M2 = float(
        stl.get("min_footprint_area_m2", base.BUILDING_MIN_FOOTPRINT_AREA_M2)
    )
    base.BUILDING_MIN_EDGE_M = float(stl.get("min_edge_m", base.BUILDING_MIN_EDGE_M))
    base.BUILDING_SIMPLIFY_M = float(stl.get("simplify_m", base.BUILDING_SIMPLIFY_M))
    height_limits = stl.get("height_limits_m")
    if height_limits and len(height_limits) == 2:
        base.BUILDING_MIN_HEIGHT_M = float(height_limits[0])
        base.BUILDING_MAX_HEIGHT_M = float(height_limits[1])


def local_to_px(x: float, y: float, radius: float) -> tuple[int, int]:
    px = PAD + int(round((x + radius) / (2.0 * radius) * (CANVAS_N - 1)))
    py = PAD + int(round((radius - y) / (2.0 * radius) * (CANVAS_N - 1)))
    return px, py


def global_ring_to_px(
    ring: list[tuple[float, float]],
    center: tuple[float, float],
    radius: float,
) -> list[tuple[int, int]]:
    return [local_to_px(x - center[0], y - center[1], radius) for x, y in ring]


def draw_inlet_arc(
    draw: ImageDraw.ImageDraw,
    radius: float,
    wind_from_deg: float,
    arc_width_deg: float,
) -> None:
    pts = []
    for bearing in np.linspace(
        wind_from_deg - 0.5 * arc_width_deg,
        wind_from_deg + 0.5 * arc_width_deg,
        360,
    ):
        theta = math.radians(float(bearing))
        pts.append(local_to_px(radius * math.sin(theta), radius * math.cos(theta), radius))
    draw.line(pts, fill=(16, 94, 205, 255), width=max(6, CANVAS_N // 180), joint="curve")


def percentile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    return ordered[int((len(ordered) - 1) * q)]


def color_gap(max_gap: float, tolerance: float, scale: float) -> tuple[int, int, int, int]:
    if max_gap <= tolerance:
        return (92, 96, 96, 75)
    t = min(1.0, max(0.0, (max_gap - tolerance) / max(scale - tolerance, 1.0e-9)))
    t = math.sqrt(t)
    # Pale yellow to deep red, so any visible warm color means an air gap.
    r0, g0, b0 = 255, 234, 153
    r1, g1, b1 = 177, 24, 38
    return (
        int(round(r0 + t * (r1 - r0))),
        int(round(g0 + t * (g1 - g0))),
        int(round(b0 + t * (b1 - b0))),
        205,
    )


def audit_buildings(
    generated: Path,
    edge_spacing_m: float,
    tolerance_m: float,
) -> tuple[list[dict], dict]:
    meta = load_metadata(generated)
    set_building_constants(meta)
    center = tuple(float(v) for v in meta["center_epsg3006"])
    radius = float(meta["radius_m"])

    land_geometry = base.load_geojson(generated / "sodermalm_osm_island_epsg3006.geojson")["features"][0]["geometry"]
    shoreline = base.polygon_boundary_rings(land_geometry)
    samples, water_level = base.terrain_samples(generated / "sodermalm_cylinder_contours.geojson", center)

    records: list[dict] = []
    skipped = 0
    for feature_index, feature in enumerate(base.load_geojson(generated / "sodermalm_buildings.geojson")["features"]):
        props = feature.get("properties", {})
        for ring_index, ring in enumerate(base.polygon_rings(feature["geometry"])):
            pts_global = base.clean_building_ring(ring)
            if len(pts_global) < 3:
                skipped += 1
                continue
            if base.BUILDING_SIMPLIFY_M > 0.0:
                pts_global = base.simplify_ring(pts_global, base.BUILDING_SIMPLIFY_M)
                pts_global = base.clean_building_ring(pts_global)
                if len(pts_global) < 3:
                    skipped += 1
                    continue

            footprint_area = abs(base.building_ring_area(pts_global))
            if footprint_area < base.BUILDING_MIN_FOOTPRINT_AREA_M2:
                skipped += 1
                continue
            if base.building_ring_min_edge(pts_global) < base.BUILDING_MIN_EDGE_M:
                skipped += 1
                continue

            cx, cy = base.ring_centroid(pts_global)
            lx, ly = cx - center[0], cy - center[1]
            if lx * lx + ly * ly > radius * radius:
                skipped += 1
                continue

            terrain_at_centroid = base.terrain_height(lx, ly, samples, water_level, shoreline, center, radius)
            if base.building_sits_on_water_level_terrain(terrain_at_centroid, water_level):
                skipped += 1
                continue

            sample_points = base.building_footprint_sample_points(pts_global, edge_spacing_m)
            terrain_values = [
                base.terrain_height(x - center[0], y - center[1], samples, water_level, shoreline, center, radius)
                for x, y in sample_points
            ]
            footprint_min_terrain = min(terrain_values)
            base_z, top_z = base.building_base_and_top(
                props, terrain_at_centroid, water_level, footprint_min_terrain
            )
            if top_z <= base_z + 0.5:
                skipped += 1
                continue

            clearances = [
                base_z - terrain_z
                for terrain_z in terrain_values
            ]
            max_gap = max(clearances)
            min_gap = min(clearances)
            records.append(
                {
                    "feature_index": feature_index,
                    "ring_index": ring_index,
                    "area_m2": footprint_area,
                    "base_z_m": base_z,
                    "top_z_m": top_z,
                    "height_m": top_z - base_z,
                    "centroid_clearance_m": base_z - terrain_at_centroid,
                    "min_clearance_m": min_gap,
                    "max_clearance_m": max_gap,
                    "has_air_gap": max_gap > tolerance_m,
                    "is_buried": min_gap < -tolerance_m,
                    "sample_count": len(sample_points),
                    "ring": pts_global,
                }
            )

    sampled_clearance_extrema = [
        value
        for record in records
        for value in (record["min_clearance_m"], record["max_clearance_m"])
    ]
    centroid_clearances = [float(record["centroid_clearance_m"]) for record in records]
    positive_gaps = [r["max_clearance_m"] for r in records if r["has_air_gap"]]
    buried_depths = [-r["min_clearance_m"] for r in records if r["is_buried"]]
    summary = {
        "generated": str(generated),
        "metadata": str(metadata_path(generated)),
        "buildings_checked": len(records),
        "buildings_skipped": skipped,
        "footprint_sample_spacing_m": edge_spacing_m,
        "tolerance_m": tolerance_m,
        "base_embed_m": base.BUILDING_BASE_EMBED_M,
        "air_gap_buildings": len(positive_gaps),
        "buried_buildings": len(buried_depths),
        "max_air_gap_m": max(positive_gaps) if positive_gaps else 0.0,
        "max_buried_depth_m": max(buried_depths) if buried_depths else 0.0,
        "sampled_clearance_m_min": min(sampled_clearance_extrema) if sampled_clearance_extrema else float("nan"),
        "sampled_clearance_m_p05": percentile(sampled_clearance_extrema, 0.05),
        "sampled_clearance_m_p50": percentile(sampled_clearance_extrema, 0.50),
        "sampled_clearance_m_p95": percentile(sampled_clearance_extrema, 0.95),
        "sampled_clearance_m_max": max(sampled_clearance_extrema) if sampled_clearance_extrema else float("nan"),
        "centroid_clearance_m_min": min(centroid_clearances) if centroid_clearances else float("nan"),
        "centroid_clearance_m_max": max(centroid_clearances) if centroid_clearances else float("nan"),
        "definition": "clearance = building_base_z - terrain_z; positive is an air gap, negative is embedded/buried",
    }
    return records, summary


def render_gap_map(
    generated: Path,
    records: list[dict],
    summary: dict,
    out_png: Path,
) -> None:
    meta = load_metadata(generated)
    center = tuple(float(v) for v in meta["center_epsg3006"])
    radius = float(meta["radius_m"])
    wind_from = float(meta.get("inflow_from_degrees", 225.0))
    arc_width = float(meta.get("inflow_arc_width_degrees", 220.0))
    scale = max(2.0, max(float(r["max_clearance_m"]) for r in records))

    canvas = Image.new("RGB", (CANVAS_N + 2 * PAD + 120, CANVAS_N + 2 * PAD), (255, 255, 255))
    disk = Image.new("RGBA", (CANVAS_N, CANVAS_N), (205, 226, 236, 255))
    disk_draw = ImageDraw.Draw(disk, "RGBA")
    disk_draw.ellipse((0, 0, CANVAS_N - 1, CANVAS_N - 1), fill=(210, 229, 218, 255))
    canvas.paste(disk.convert("RGB"), (PAD, PAD))
    draw = ImageDraw.Draw(canvas, "RGBA")

    land_path = generated / "sodermalm_osm_island_epsg3006.geojson"
    if land_path.exists():
        land = base.load_geojson(land_path)["features"][0]["geometry"]
        for ring in base.polygon_boundary_rings(land):
            draw.line(global_ring_to_px(ring, center, radius), fill=(15, 25, 28, 220), width=max(2, CANVAS_N // 520), joint="curve")

    for record in sorted(records, key=lambda item: item["max_clearance_m"]):
        pts = global_ring_to_px(record["ring"], center, radius)
        if len(pts) < 3:
            continue
        fill = color_gap(float(record["max_clearance_m"]), float(summary["tolerance_m"]), scale)
        outline = (38, 64, 74, 95)
        if record["is_buried"]:
            outline = (32, 73, 171, 190)
        draw.polygon(pts, fill=fill, outline=outline)

    draw_inlet_arc(draw, radius, wind_from, arc_width)

    bar_x = PAD + CANVAS_N + 32
    bar_y = PAD + CANVAS_N // 8
    bar_h = 3 * CANVAS_N // 4
    bar_w = 32
    for iy in range(bar_h):
        gap = scale * (1.0 - iy / max(bar_h - 1, 1))
        draw.rectangle(
            (bar_x, bar_y + iy, bar_x + bar_w, bar_y + iy),
            fill=color_gap(gap, float(summary["tolerance_m"]), scale),
        )
    draw.rectangle((bar_x, bar_y, bar_x + bar_w, bar_y + bar_h), outline=(20, 20, 20, 255), width=1)
    for frac, label in ((0.0, f"{scale:.1f}"), (0.5, f"{0.5 * scale:.1f}"), (1.0, "0")):
        yy = bar_y + int(round(frac * bar_h))
        draw.line((bar_x + bar_w, yy, bar_x + bar_w + 8, yy), fill=(20, 20, 20, 255), width=1)
        draw.text((bar_x + bar_w + 12, yy - 7), label, fill=(20, 20, 20, 255))
    draw.text((bar_x, bar_y - 24), "gap m", fill=(20, 20, 20, 255))

    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_png)


def serializable_records(records: list[dict]) -> list[dict]:
    cleaned = []
    for record in records:
        item = dict(record)
        item.pop("ring", None)
        cleaned.append(item)
    return cleaned


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("generated", nargs="?", type=Path, default=DEFAULT_GENERATED)
    parser.add_argument("--edge-spacing-m", type=float, default=5.0)
    parser.add_argument("--tolerance-m", type=float, default=0.05)
    parser.add_argument("--out-dir", type=Path)
    args = parser.parse_args()

    generated = args.generated.resolve()
    if args.out_dir is None:
        out_dir = generated.parent / f"figures_{generated.name.removeprefix('generated_')}" / "building_ground_contact"
    else:
        out_dir = args.out_dir.resolve()

    records, summary = audit_buildings(generated, args.edge_spacing_m, args.tolerance_m)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "building_ground_contact_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    worst = sorted(serializable_records(records), key=lambda item: item["max_clearance_m"], reverse=True)[:200]
    (out_dir / "building_ground_contact_worst_gaps.json").write_text(json.dumps(worst, indent=2), encoding="utf-8")
    render_gap_map(generated, records, summary, out_dir / "building_ground_contact_gaps_topdown.png")

    print(json.dumps(summary, indent=2))
    print(out_dir / "building_ground_contact_gaps_topdown.png")


if __name__ == "__main__":
    main()
