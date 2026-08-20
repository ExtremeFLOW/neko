#!/usr/bin/env python3
"""Fetch and clean the OSM island polygon for the Södermalm mesh.

The Nominatim `place=island` polygon is too generous for this CFD cutout: it
can include land around Hammarby sjö, e.g. Lumaparken/Hammarby Sjöstad.  We
therefore use OSM water polygons as cutters and keep only the connected
component containing central Södermalm.
"""

from __future__ import annotations

import json
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path

from osgeo import ogr, osr


HERE = Path(__file__).resolve().parent
OUT = HERE / "generated"
OGR2OGR = "/opt/homebrew/bin/ogr2ogr"
URL = "https://nominatim.openstreetmap.org/search?q=S%C3%B6dermalm%2C%20Stockholm&format=json&polygon_geojson=1&limit=5"
WATER_QUERIES = [
    "Hammarby sjö, Stockholm",
    "Årstaviken, Stockholm",
    "Riddarfjärden, Stockholm",
]
CENTRAL_SODERMALM_WGS84 = (18.063, 59.314)
WATER_CUT_BUFFER_M = 18.0


def fetch_nominatim(query: str) -> list[dict]:
    url = "https://nominatim.openstreetmap.org/search?" + urllib.parse.urlencode(
        {"q": query, "format": "json", "polygon_geojson": "1", "limit": "5"}
    )
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Codex local geometry test for Neko Södermalm mesh"},
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.loads(response.read().decode("utf-8"))


def epsg_transform(src_epsg: int, dst_epsg: int) -> osr.CoordinateTransformation:
    src = osr.SpatialReference()
    src.ImportFromEPSG(src_epsg)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(dst_epsg)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return osr.CoordinateTransformation(src, dst)


def geojson_to_epsg3006(geometry: dict) -> ogr.Geometry:
    geom = ogr.CreateGeometryFromJson(json.dumps(geometry))
    if geom is None:
        raise SystemExit("Could not parse OSM geometry.")
    geom.Transform(epsg_transform(4326, 3006))
    return geom.MakeValid()


def iter_polygons(geometry: ogr.Geometry) -> list[ogr.Geometry]:
    if geometry.GetGeometryName() == "POLYGON":
        return [geometry.Clone()]
    if geometry.GetGeometryName() in {"MULTIPOLYGON", "GEOMETRYCOLLECTION"}:
        parts: list[ogr.Geometry] = []
        for i in range(geometry.GetGeometryCount()):
            parts.extend(iter_polygons(geometry.GetGeometryRef(i)))
        return parts
    return []


def write_geojson(path: Path, geometry: ogr.Geometry, properties: dict, epsg: int = 3006) -> None:
    if path.exists():
        path.unlink()
    driver = ogr.GetDriverByName("GeoJSON")
    ds = driver.CreateDataSource(str(path))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    layer = ds.CreateLayer(path.stem, srs=srs, geom_type=geometry.GetGeometryType())
    for key, value in properties.items():
        field = ogr.FieldDefn(key, ogr.OFTString)
        layer.CreateField(field)
    feature = ogr.Feature(layer.GetLayerDefn())
    for key, value in properties.items():
        feature.SetField(key, str(value))
    feature.SetGeometry(geometry)
    layer.CreateFeature(feature)
    feature = None
    ds = None


def point_epsg3006(lon: float, lat: float) -> ogr.Geometry:
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(lon, lat)
    point.Transform(epsg_transform(4326, 3006))
    return point


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wgs = OUT / "sodermalm_osm_island_wgs84.geojson"
    epsg3006_raw = OUT / "sodermalm_osm_island_raw_epsg3006.geojson"
    epsg3006 = OUT / "sodermalm_osm_island_epsg3006.geojson"
    water_cut = OUT / "sodermalm_osm_water_cut_epsg3006.geojson"

    req = urllib.request.Request(URL, headers={"User-Agent": "Codex local geometry test for Neko Södermalm mesh"})
    with urllib.request.urlopen(req, timeout=30) as response:
        data = json.loads(response.read().decode("utf-8"))

    island = None
    for item in data:
        if item.get("class") == "place" and item.get("type") == "island" and item.get("geojson", {}).get("type") in {"Polygon", "MultiPolygon"}:
            island = item
            break
    if island is None:
        raise SystemExit("Could not find OSM place=island polygon for Södermalm.")

    feature_collection = {
        "type": "FeatureCollection",
        "name": wgs.stem,
        "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "source": "OpenStreetMap Nominatim",
                    "osm_type": island.get("osm_type"),
                    "osm_id": island.get("osm_id"),
                    "class": island.get("class"),
                    "type": island.get("type"),
                    "display_name": island.get("display_name"),
                    "licence": island.get("licence"),
                },
                "geometry": island["geojson"],
            }
        ],
    }
    wgs.write_text(json.dumps(feature_collection, indent=2), encoding="utf-8")
    subprocess.run([OGR2OGR, "-f", "GeoJSON", "-t_srs", "EPSG:3006", str(epsg3006_raw), str(wgs)], check=True)

    island_geom = geojson_to_epsg3006(island["geojson"])
    water_union = None
    water_props = []
    for query in WATER_QUERIES:
        match = None
        for item in fetch_nominatim(query):
            if item.get("geojson", {}).get("type") in {"Polygon", "MultiPolygon"}:
                match = item
                break
        if match is None:
            raise SystemExit(f"Could not find OSM water polygon for {query}.")
        geom = geojson_to_epsg3006(match["geojson"])
        water_union = geom.Clone() if water_union is None else water_union.Union(geom)
        water_props.append(f'{match.get("display_name")} [{match.get("osm_type")} {match.get("osm_id")}]')

    assert water_union is not None
    water_union = water_union.MakeValid().Buffer(WATER_CUT_BUFFER_M).MakeValid()
    cut = island_geom.Difference(water_union).MakeValid()
    center_point = point_epsg3006(*CENTRAL_SODERMALM_WGS84)
    parts = iter_polygons(cut)
    if not parts:
        raise SystemExit("Water-cut Södermalm mask produced no polygons.")
    selected = None
    for part in parts:
        if part.Contains(center_point):
            selected = part
            break
    if selected is None:
        selected = max(parts, key=lambda geom: geom.GetArea())

    props = dict(feature_collection["features"][0]["properties"])
    props.update(
        {
            "cleaning": "raw OSM polygon minus buffered OSM water polygons; selected central Södermalm component",
            "water_cut_buffer_m": WATER_CUT_BUFFER_M,
            "water_cut_sources": " | ".join(water_props),
            "raw_area_m2": f"{island_geom.GetArea():.1f}",
            "clean_area_m2": f"{selected.GetArea():.1f}",
        }
    )
    write_geojson(water_cut, water_union, {"sources": " | ".join(water_props), "buffer_m": WATER_CUT_BUFFER_M})
    write_geojson(epsg3006, selected, props)
    print(json.dumps(props, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
