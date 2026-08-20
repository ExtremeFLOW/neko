#!/usr/bin/env python3
"""Derive a Södermalm land mask from local contour and building geometry."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from osgeo import ogr, osr


HERE = Path(__file__).resolve().parent
INPUT_GPKG = HERE / "input" / "inner_stockholm_15km_layers_clipped.gpkg"
OUT = HERE / "generated"

CONTOUR_BUFFER_M = 35.0
BUILDING_BUFFER_M = 12.0
CLOSING_BUFFER_M = 45.0
SIMPLIFY_M = 12.0

SODERMALM_WGS84_HINT = [
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


def transform_hint_to_epsg3006() -> ogr.Geometry:
    src = osr.SpatialReference()
    src.ImportFromEPSG(4326)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(3006)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transform = osr.CoordinateTransformation(src, dst)

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for lon, lat in SODERMALM_WGS84_HINT:
        x, y, _ = transform.TransformPoint(lon, lat)
        ring.AddPoint(x, y)
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


def bbox_geometry(bbox: tuple[float, float, float, float]) -> ogr.Geometry:
    xmin, ymin, xmax, ymax = bbox
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in ((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)):
        ring.AddPoint(x, y)
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


def iter_layer_geoms(ds, layer_name: str, bbox: tuple[float, float, float, float]):
    clip = bbox_geometry(bbox)
    layer = ds.GetLayerByName(layer_name)
    layer.SetSpatialFilterRect(*bbox)
    layer.ResetReading()
    for feature in layer:
        geom = feature.GetGeometryRef()
        if geom:
            clipped = geom.Intersection(clip)
            if clipped and not clipped.IsEmpty():
                yield clipped
    layer.SetSpatialFilter(None)


def add_buffered(collection: ogr.Geometry, geom: ogr.Geometry, distance: float) -> None:
    if geom.IsEmpty():
        return
    buffered = geom.Buffer(distance, 4)
    if buffered and not buffered.IsEmpty():
        for part in parts(buffered):
            collection.AddGeometry(part)


def parts(geom: ogr.Geometry) -> list[ogr.Geometry]:
    gtype = geom.GetGeometryType()
    flat = ogr.GT_Flatten(gtype)
    if flat == ogr.wkbPolygon:
        return [geom]
    if flat == ogr.wkbMultiPolygon or flat == ogr.wkbGeometryCollection:
        return [geom.GetGeometryRef(i).Clone() for i in range(geom.GetGeometryCount()) if ogr.GT_Flatten(geom.GetGeometryRef(i).GetGeometryType()) == ogr.wkbPolygon]
    return []


def geom_to_geojson_polygon(geom: ogr.Geometry) -> dict:
    return json.loads(geom.ExportToJson())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ds = ogr.Open(str(INPUT_GPKG))
    if ds is None:
        raise SystemExit(f"Could not open {INPUT_GPKG}")

    hint = transform_hint_to_epsg3006()
    xmin, xmax, ymin, ymax = hint.GetEnvelope()
    # Expand enough to include the water/shoreline immediately around the island,
    # but not so much that nearby islands dominate the connected component.
    pad = 420.0
    bbox = (xmin - pad, ymin - pad, xmax + pad, ymax + pad)
    seed = hint.Centroid()

    collection = ogr.Geometry(ogr.wkbMultiPolygon)
    n_contours = 0
    for geom in iter_layer_geoms(ds, "terrain_contours", bbox):
        add_buffered(collection, geom, CONTOUR_BUFFER_M)
        n_contours += 1

    n_buildings = 0
    for geom in iter_layer_geoms(ds, "buildings", bbox):
        add_buffered(collection, geom, BUILDING_BUFFER_M)
        n_buildings += 1

    union = collection.UnionCascaded()
    if union is None or union.IsEmpty():
        raise SystemExit("Could not derive a land union from contours/buildings.")

    # Morphological close/open: fill block-scale holes in the land component
    # without spanning broad waterways.
    cleaned = union.Buffer(CLOSING_BUFFER_M, 4).Buffer(-CLOSING_BUFFER_M, 4).MakeValid()
    candidates = []
    for part in parts(cleaned):
        if part.IsEmpty():
            continue
        contains_seed = part.Contains(seed)
        overlap = part.Intersection(hint).Area()
        candidates.append((contains_seed, overlap, part.Area(), part))
    if not candidates:
        raise SystemExit("Derived land union contains no polygon candidates.")

    candidates.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
    selected = candidates[0][3].SimplifyPreserveTopology(SIMPLIFY_M).MakeValid()
    if ogr.GT_Flatten(selected.GetGeometryType()) != ogr.wkbPolygon:
        selected_parts = parts(selected)
        selected_parts.sort(key=lambda g: g.Area(), reverse=True)
        selected = selected_parts[0]

    out_path = OUT / "sodermalm_land_mask_derived.geojson"
    feature_collection = {
        "type": "FeatureCollection",
        "name": out_path.stem,
        "crs": {"type": "name", "properties": {"name": "EPSG:3006"}},
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "method": "contour_building_buffer_component",
                    "contour_buffer_m": CONTOUR_BUFFER_M,
                    "building_buffer_m": BUILDING_BUFFER_M,
                    "closing_buffer_m": CLOSING_BUFFER_M,
                    "simplify_m": SIMPLIFY_M,
                    "source_contours": n_contours,
                    "source_buildings": n_buildings,
                    "area_m2": selected.Area(),
                    "perimeter_m": selected.Boundary().Length(),
                },
                "geometry": geom_to_geojson_polygon(selected),
            }
        ],
    }
    out_path.write_text(json.dumps(feature_collection, indent=2), encoding="utf-8")
    print(json.dumps(feature_collection["features"][0]["properties"], indent=2))


if __name__ == "__main__":
    main()
