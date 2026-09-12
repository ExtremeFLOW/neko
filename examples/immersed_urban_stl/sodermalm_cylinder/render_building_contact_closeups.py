#!/usr/bin/env python3
"""Render close-up 3D checks of building bases relative to terrain."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

import build_sodermalm_cylinder as geom
from check_building_ground_contact import audit_buildings, load_metadata, set_building_constants


HERE = Path(__file__).resolve().parent
DEFAULT_GENERATED = HERE / "generated_p7_xy35_inlet220_nz10"
W, H = 1400, 1000


def terrain_color(z: float, zmin: float, zmax: float, alpha: int = 120) -> tuple[int, int, int, int]:
    t = 0.0 if zmax <= zmin else (z - zmin) / (zmax - zmin)
    t = max(0.0, min(1.0, t))
    lo = np.array([198, 218, 205], dtype=float)
    hi = np.array([148, 180, 151], dtype=float)
    rgb = lo + t * (hi - lo)
    return (int(rgb[0]), int(rgb[1]), int(rgb[2]), alpha)


def shade(rgb: tuple[int, int, int], factor: float, alpha: int = 230) -> tuple[int, int, int, int]:
    return (
        max(0, min(255, int(rgb[0] * factor))),
        max(0, min(255, int(rgb[1] * factor))),
        max(0, min(255, int(rgb[2] * factor))),
        alpha,
    )


class Camera:
    def __init__(self, bounds: tuple[float, float, float, float, float, float]) -> None:
        xmin, xmax, ymin, ymax, zmin, zmax = bounds
        self.cx = 0.5 * (xmin + xmax)
        self.cy = 0.5 * (ymin + ymax)
        self.cz = 0.5 * (zmin + zmax)
        self.yaw = math.radians(-38.0)
        self.pitch = math.radians(28.0)
        xy_span = max(xmax - xmin, ymax - ymin, 1.0)
        z_span = max(zmax - zmin, 1.0)
        self.z_scale = min(2.0, max(1.0, xy_span / (5.0 * z_span)))
        self.scale = 0.70 * min(W, H) / max(xy_span, 2.2 * z_span * self.z_scale)

    def project(self, p: tuple[float, float, float]) -> tuple[float, float, float]:
        x = p[0] - self.cx
        y = p[1] - self.cy
        z = (p[2] - self.cz) * self.z_scale
        ca, sa = math.cos(self.yaw), math.sin(self.yaw)
        x1 = ca * x - sa * y
        y1 = sa * x + ca * y
        cp, sp = math.cos(self.pitch), math.sin(self.pitch)
        y2 = cp * y1 - sp * z
        z2 = sp * y1 + cp * z
        return (W * 0.5 + self.scale * x1, H * 0.57 - self.scale * y2, z2)


def polygon_area(points: list[tuple[float, float]]) -> float:
    return 0.5 * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(points, points[1:] + points[:1]))


def render_record(generated: Path, record: dict, out: Path, title: str) -> None:
    meta = load_metadata(generated)
    set_building_constants(meta)
    center = tuple(float(v) for v in meta["center_epsg3006"])
    radius = float(meta["radius_m"])
    land_geometry = geom.load_geojson(generated / "sodermalm_osm_island_epsg3006.geojson")["features"][0]["geometry"]
    shoreline = geom.polygon_boundary_rings(land_geometry)
    samples, water_level = geom.terrain_samples(generated / "sodermalm_cylinder_contours.geojson", center)

    ring_global = record["ring"]
    ring = [(x - center[0], y - center[1]) for x, y in ring_global]
    if polygon_area(ring) < 0.0:
        ring.reverse()

    base_z = float(record["base_z_m"])
    top_z = float(record["top_z_m"])
    xs = [p[0] for p in ring]
    ys = [p[1] for p in ring]
    span = max(max(xs) - min(xs), max(ys) - min(ys), 40.0)
    margin = max(25.0, 0.8 * span)
    xmin, xmax = min(xs) - margin, max(xs) + margin
    ymin, ymax = min(ys) - margin, max(ys) + margin

    nx = ny = 74
    gx = np.linspace(xmin, xmax, nx)
    gy = np.linspace(ymin, ymax, ny)
    z = np.zeros((ny, nx), dtype=float)
    for j, y in enumerate(gy):
        for i, x in enumerate(gx):
            z[j, i] = geom.terrain_height(x, y, samples, water_level, shoreline, center, radius)

    zmin = min(float(np.min(z)), base_z)
    zmax = max(float(np.max(z)), top_z)
    cam = Camera((xmin, xmax, ymin, ymax, zmin, zmax))

    img = Image.new("RGBA", (W, H), (238, 244, 247, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    items: list[tuple[float, list[tuple[float, float]], tuple[int, int, int, int], tuple[int, int, int, int] | None]] = []

    terrain_zmin = float(np.min(z))
    terrain_zmax = float(np.max(z))
    for j in range(ny - 1):
        for i in range(nx - 1):
            verts = [
                (float(gx[i]), float(gy[j]), float(z[j, i])),
                (float(gx[i + 1]), float(gy[j]), float(z[j, i + 1])),
                (float(gx[i + 1]), float(gy[j + 1]), float(z[j + 1, i + 1])),
                (float(gx[i]), float(gy[j + 1]), float(z[j + 1, i])),
            ]
            pts3 = [cam.project(v) for v in verts]
            depth = sum(p[2] for p in pts3) / 4.0
            zc = sum(v[2] for v in verts) / 4.0
            items.append((depth, [(p[0], p[1]) for p in pts3], terrain_color(zc, terrain_zmin, terrain_zmax), None))

    roof = [(x, y, top_z) for x, y in ring]
    bottom = [(x, y, base_z) for x, y in ring]
    wall_rgb = (142, 146, 148)
    for i in range(len(ring)):
        j = (i + 1) % len(ring)
        verts = [bottom[i], bottom[j], roof[j], roof[i]]
        pts3 = [cam.project(v) for v in verts]
        depth = sum(p[2] for p in pts3) / 4.0
        dx = ring[j][0] - ring[i][0]
        dy = ring[j][1] - ring[i][1]
        factor = 0.72 + 0.24 * abs(math.sin(math.atan2(dy, dx) - cam.yaw))
        items.append((depth, [(p[0], p[1]) for p in pts3], shade(wall_rgb, factor, 245), (76, 82, 85, 170)))

    roof_pts3 = [cam.project(v) for v in roof]
    bottom_pts3 = [cam.project(v) for v in bottom]
    items.append((sum(p[2] for p in bottom_pts3) / len(bottom_pts3), [(p[0], p[1]) for p in bottom_pts3], (96, 100, 104, 165), (55, 60, 64, 190)))
    items.append((sum(p[2] for p in roof_pts3) / len(roof_pts3), [(p[0], p[1]) for p in roof_pts3], (184, 187, 188, 255), (92, 96, 98, 170)))

    for _, pts, fill, outline in sorted(items, key=lambda item: item[0]):
        draw.polygon(pts, fill=fill)
        if outline is not None:
            draw.line(pts + [pts[0]], fill=outline, width=1)

    base_outline = [cam.project(p) for p in bottom]
    terrain_outline = [
        cam.project(
            (
                x,
                y,
                geom.terrain_height(x, y, samples, water_level, shoreline, center, radius),
            )
        )
        for x, y in ring
    ]
    draw.line(
        [(p[0], p[1]) for p in base_outline] + [(base_outline[0][0], base_outline[0][1])],
        fill=(180, 35, 35, 235),
        width=3,
        joint="curve",
    )
    draw.line(
        [(p[0], p[1]) for p in terrain_outline] + [(terrain_outline[0][0], terrain_outline[0][1])],
        fill=(0, 125, 150, 240),
        width=3,
        joint="curve",
    )

    # Draw a small vertical ruler at the sample point with the deepest burial.
    sample_points = geom.building_footprint_sample_points(ring_global, float(meta["building_stl"]["contact_sample_spacing_m"]))
    terrain_values = [
        geom.terrain_height(x - center[0], y - center[1], samples, water_level, shoreline, center, radius)
        for x, y in sample_points
    ]
    k = int(np.argmax(terrain_values))
    sx, sy = sample_points[k][0] - center[0], sample_points[k][1] - center[1]
    st = float(terrain_values[k])
    p0 = cam.project((sx, sy, base_z))
    p1 = cam.project((sx, sy, st))
    draw.line((p0[0], p0[1], p1[0], p1[1]), fill=(190, 40, 40, 230), width=4)
    r = 5
    draw.ellipse((p1[0] - r, p1[1] - r, p1[0] + r, p1[1] + r), fill=(190, 40, 40, 235))

    draw.text((34, 28), title, fill=(20, 25, 28, 255))
    draw.text(
        (34, 56),
        f"base embed range: {-record['max_clearance_m']:.2f} to {-record['min_clearance_m']:.2f} m",
        fill=(20, 25, 28, 255),
    )
    draw.text((34, 84), "red marker: deepest local burial", fill=(120, 35, 35, 255))
    draw.text((34, 112), "red line: building base; teal line: terrain at wall", fill=(20, 70, 80, 255))

    out.parent.mkdir(parents=True, exist_ok=True)
    img.convert("RGB").save(out)
    print(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("generated", nargs="?", type=Path, default=DEFAULT_GENERATED)
    parser.add_argument("--out-dir", type=Path)
    args = parser.parse_args()

    generated = args.generated.resolve()
    records, _ = audit_buildings(generated, edge_spacing_m=5.0, tolerance_m=0.05)
    out_dir = args.out_dir or (
        generated.parent / f"figures_{generated.name.removeprefix('generated_')}" / "building_ground_contact"
    )
    most = min(records, key=lambda record: record["min_clearance_m"])
    visible = [record for record in records if record["area_m2"] >= 200.0]
    least = max(visible or records, key=lambda record: record["min_clearance_m"])

    render_record(generated, most, out_dir / "building_contact_most_buried_3d.png", "Most buried building")
    render_record(generated, least, out_dir / "building_contact_least_buried_3d.png", "Least buried building")


if __name__ == "__main__":
    main()
