#!/usr/bin/env python3
"""Expand the prepared Södermalm cutout with the Gamla Stan island group."""

from __future__ import annotations

import json
import math
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path


HERE = Path(__file__).resolve().parent
INPUT = HERE / "input"
PREPARED = INPUT / "prepared"
WORK = INPUT / "prepared_work"
INPUT_GPKG = INPUT / "inner_stockholm_15km_layers_clipped.gpkg"
OGR2OGR = "/opt/homebrew/bin/ogr2ogr"

LAND_MASK = PREPARED / "sodermalm_osm_island_epsg3006.geojson"
SODERMALM_ONLY_MASK = PREPARED / "sodermalm_only_epsg3006.geojson"
BUILDINGS = PREPARED / "sodermalm_buildings.geojson"
CONTOURS = PREPARED / "sodermalm_cylinder_contours.geojson"

UPSTREAM_ISLAND_QUERIES = (
    "Stadsholmen, Stockholm",
    "Riddarholmen, Stockholm",
)

BUFFER_M = 220.0
UPSTREAM_CLOSE_M = 0.0


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def fetch_nominatim(query: str) -> list[dict]:
    url = "https://nominatim.openstreetmap.org/search?" + urllib.parse.urlencode(
        {"q": query, "format": "json", "polygon_geojson": "1", "limit": "5"}
    )
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Codex local geometry test for Neko Stockholm mesh"},
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.loads(response.read().decode("utf-8"))


def pick_polygon(query: str) -> dict:
    for item in fetch_nominatim(query):
        geometry = item.get("geojson", {})
        if geometry.get("type") in {"Polygon", "MultiPolygon"}:
            return {
                "type": "Feature",
                "properties": {
                    "query": query,
                    "display_name": item.get("display_name", ""),
                    "osm_type": item.get("osm_type", ""),
                    "osm_id": item.get("osm_id", ""),
                    "class": item.get("class", ""),
                    "type": item.get("type", ""),
                },
                "geometry": geometry,
            }
    raise SystemExit(f"Could not find a polygon result for {query}")


def write_geojson(path: Path, features: list[dict], crs: str) -> None:
    data = {
        "type": "FeatureCollection",
        "name": path.stem,
        "crs": {"type": "name", "properties": {"name": crs}},
        "features": features,
    }
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def iter_xy(geometry: dict):
    gtype = geometry.get("type")
    coords = geometry.get("coordinates", [])
    if gtype == "Polygon":
        for ring in coords:
            for x, y, *_ in ring:
                yield float(x), float(y)
    elif gtype == "MultiPolygon":
        for polygon in coords:
            for ring in polygon:
                for x, y, *_ in ring:
                    yield float(x), float(y)


def polygon_area(ring: list[list[float]]) -> float:
    area = 0.0
    for p0, p1 in zip(ring, ring[1:] + ring[:1]):
        area += float(p0[0]) * float(p1[1]) - float(p1[0]) * float(p0[1])
    return 0.5 * area


def extract_largest_polygon(source: Path, target: Path) -> None:
    data = json.loads(source.read_text(encoding="utf-8"))
    polygons: list[list] = []
    for feature in data["features"]:
        geom = feature["geometry"]
        if geom["type"] == "Polygon":
            polygons.append(geom["coordinates"])
        elif geom["type"] == "MultiPolygon":
            polygons.extend(geom["coordinates"])
    if not polygons:
        raise SystemExit(f"No polygon found in {source}")
    selected = max(polygons, key=lambda poly: abs(polygon_area(poly[0])))
    write_geojson(
        target,
        [
            {
                "type": "Feature",
                "properties": {"source": source.name, "selection": "largest_polygon_sodermalm_only"},
                "geometry": {"type": "Polygon", "coordinates": selected},
            }
        ],
        "EPSG:3006",
    )


def land_bbox(path: Path) -> tuple[float, float, float, float]:
    data = json.loads(path.read_text(encoding="utf-8"))
    points = [xy for feature in data["features"] for xy in iter_xy(feature["geometry"])]
    if not points:
        raise SystemExit(f"No polygon coordinates found in {path}")
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    center = ((min(xs) + max(xs)) * 0.5, (min(ys) + max(ys)) * 0.5)
    radius = max(math.hypot(x - center[0], y - center[1]) for x, y in points) + BUFFER_M
    return center[0] - radius, center[1] - radius, center[0] + radius, center[1] + radius


def main() -> None:
    PREPARED.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    if not SODERMALM_ONLY_MASK.exists():
        if not LAND_MASK.exists():
            raise SystemExit(f"Missing base Södermalm land mask: {LAND_MASK}")
        extract_largest_polygon(LAND_MASK, SODERMALM_ONLY_MASK)
    if not INPUT_GPKG.exists():
        raise SystemExit(f"Missing input GeoPackage: {INPUT_GPKG}")

    upstream_wgs = WORK / "gamla_stan_islands_wgs84.geojson"
    upstream_epsg3006 = WORK / "gamla_stan_islands_epsg3006.geojson"
    expanded_gpkg = WORK / "expanded_land.gpkg"
    expanded_land = WORK / "expanded_land_epsg3006.geojson"

    write_geojson(upstream_wgs, [pick_polygon(query) for query in UPSTREAM_ISLAND_QUERIES], "EPSG:4326")
    for path in (upstream_epsg3006, expanded_land):
        if path.exists():
            path.unlink()
    run([OGR2OGR, "-overwrite", "-f", "GeoJSON", "-t_srs", "EPSG:3006", str(upstream_epsg3006), str(upstream_wgs)])

    if expanded_gpkg.exists():
        expanded_gpkg.unlink()
    run([OGR2OGR, "-overwrite", "-f", "GPKG", str(expanded_gpkg), str(SODERMALM_ONLY_MASK), "-nln", "sodermalm", "-nlt", "PROMOTE_TO_MULTI"])
    run([OGR2OGR, "-update", str(expanded_gpkg), str(upstream_epsg3006), "-nln", "upstream", "-nlt", "PROMOTE_TO_MULTI"])
    run([
        OGR2OGR,
        "-overwrite",
        "-f",
        "GeoJSON",
        str(expanded_land),
        str(expanded_gpkg),
        "-dialect",
        "sqlite",
        "-sql",
        (
            "SELECT ST_Union(geom) AS geom FROM ("
            "SELECT geom FROM sodermalm "
            "UNION ALL "
            "SELECT ST_Union(geom) AS geom FROM upstream"
            ")"
        ),
        "-nln",
        LAND_MASK.stem,
    ])
    expanded_land.replace(LAND_MASK)

    xmin, ymin, xmax, ymax = land_bbox(LAND_MASK)
    for path in (BUILDINGS, CONTOURS):
        if path.exists():
            path.unlink()
    run([OGR2OGR, "-f", "GeoJSON", str(BUILDINGS), str(INPUT_GPKG), "buildings", "-clipsrc", str(LAND_MASK)])
    run([
        OGR2OGR,
        "-f",
        "GeoJSON",
        str(CONTOURS),
        str(INPUT_GPKG),
        "terrain_contours",
        "-spat",
        str(xmin),
        str(ymin),
        str(xmax),
        str(ymax),
    ])

    print(
        json.dumps(
            {
                "land_mask": str(LAND_MASK),
                "buildings": str(BUILDINGS),
                "contours": str(CONTOURS),
                "added_islands": UPSTREAM_ISLAND_QUERIES,
                "upstream_close_m": UPSTREAM_CLOSE_M,
                "upstream_shape_note": (
                    "Uses the Stadsholmen islet polygon instead of the broader "
                    "Gamla stan administrative polygon; no manual southeast cut "
                    "and no buffer-based closing."
                ),
                "bbox_epsg3006": [xmin, ymin, xmax, ymax],
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
