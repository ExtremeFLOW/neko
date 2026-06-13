#!/usr/bin/env python3
"""Create a Brinkman mask directly from GIS building footprints.

Neko's boundary_mesh signed-distance path expects a well-behaved closed surface.
The urban building STL is a collection of many disconnected building solids,
which can produce false solid regions. This script avoids that ambiguity by
marking GLL points that are inside the GIS footprint and between the building
base and roof elevations.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pymech.neksuite import readnek, writenek
from shapely import box, contains_xy
from shapely.geometry import MultiPolygon, Polygon, shape
from shapely.strtree import STRtree


HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parents[1]
GENERATED = EXAMPLE / "geometry_pipeline" / "generated"
TEMPLATE = EXAMPLE / "fine_3x" / "run_mask" / "fields" / "brinkman0.f00000"
OUTPUT_DIR = EXAMPLE / "fine_3x" / "upload"
OUTPUT_NAME = "urban_brinkman_mask"


def local_polygon(poly: Polygon, origin_x: float, origin_y: float) -> Polygon:
    exterior = [(p[0] - origin_x, p[1] - origin_y) for p in poly.exterior.coords]
    holes = [
        [(p[0] - origin_x, p[1] - origin_y) for p in ring.coords]
        for ring in poly.interiors
    ]
    result = Polygon(exterior, holes)
    if not result.is_valid:
        result = result.buffer(0)
    return result


def building_parts() -> tuple[list[Polygon], np.ndarray, np.ndarray]:
    metadata = json.loads((GENERATED / "metadata.json").read_text())
    origin_x = metadata["local_origin_epsg3006"]["x"]
    origin_y = metadata["local_origin_epsg3006"]["y"]
    geojson = json.loads((GENERATED / "buildings_epsg3006.geojson").read_text())

    polygons: list[Polygon] = []
    bases: list[float] = []
    tops: list[float] = []

    for feature in geojson["features"]:
        geom = shape(feature["geometry"])
        props = feature.get("properties", {})
        base = float(props.get("MARK_Z") or props.get("TAKMIN_Z") or 0.0)
        top = float(props.get("TAKMAX_Z") or props.get("TAK_Z") or 0.0)
        if top <= base:
            top = base + float(props.get("BYGG_H") or 0.0)
        if top <= base:
            continue

        parts = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
        for part in parts:
            if not isinstance(part, Polygon):
                continue
            polygon = local_polygon(part, origin_x, origin_y)
            if polygon.is_empty:
                continue
            polygons.append(polygon)
            bases.append(base)
            tops.append(top)

    return polygons, np.asarray(bases), np.asarray(tops)


def main() -> None:
    polygons, bases, tops = building_parts()
    tree = STRtree(polygons)

    data = readnek(str(TEMPLATE))
    solid_points = 0
    total_points = 0

    for elem in data.elem:
        x = elem.pos[0].ravel()
        y = elem.pos[1].ravel()
        z = elem.pos[2].ravel()
        mask = np.zeros(x.shape, dtype=float)

        candidates = tree.query(
            box(float(x.min()), float(y.min()), float(x.max()), float(y.max()))
        )
        for idx in np.asarray(candidates, dtype=int):
            inside = contains_xy(polygons[idx], x, y)
            if not np.any(inside):
                continue
            mask = np.maximum(mask, inside & (z >= bases[idx]) & (z <= tops[idx]))

        scalars = np.asarray(elem.scal)
        scalars[0].ravel()[:] = mask
        if scalars.shape[0] > 1:
            scalars[1].ravel()[:] = 1000.0 * mask

        solid_points += int(mask.sum())
        total_points += mask.size

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    writenek(str(OUTPUT_DIR / f"{OUTPUT_NAME}0.f00000"), data)
    (OUTPUT_DIR / f"{OUTPUT_NAME}0.nek5000").write_text(
        f"filetemplate:         {OUTPUT_NAME}%01d.f%05d\n"
        "firsttimestep:     0\n"
        "numtimesteps:     1\n"
    )
    print(f"Solid fraction: {solid_points / total_points:.6f}")
    print(f"Wrote {OUTPUT_DIR / f'{OUTPUT_NAME}0.f00000'}")


if __name__ == "__main__":
    main()
