#!/usr/bin/env python3
"""Build a 2x finer horizontal Södermalm cylinder mesh.

This keeps the validated coarse-run vertical setup but halves the x-y target
element sizes to give roughly 55 m effective p2 spacing over land.
"""

from __future__ import annotations

import json
import math
import os
import shutil
from pathlib import Path

import numpy as np

import build_sodermalm_cylinder as base


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "input" / "prepared"
OUT = HERE / "generated_2x_xy"
FIG = HERE / "figures_2x_xy"
TAG = "2x_xy"

PREPARED_INPUTS = (
    "sodermalm_osm_island_epsg3006.geojson",
    "sodermalm_osm_water_cut_epsg3006.geojson",
    "sodermalm_buildings.geojson",
    "sodermalm_cylinder_contours.geojson",
)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    skip_figures = os.environ.get("SODERMALM_SKIP_FIGURES", "0") == "1"

    base.OUT = OUT
    base.FIG = FIG
    base.NZ = 5
    base.DOMAIN_HEIGHT_ABOVE_LOWEST_M = 350.0
    base.VERTICAL_STRETCH = 2.0
    base.LAND_MESH_SIZE_M = 110.0
    base.WATER_MESH_SIZE_M = 225.0
    base.SHORELINE_SIMPLIFY_M = 35.0
    base.SHORE_BLEND_M = 180.0

    missing = [name for name in PREPARED_INPUTS if not (SOURCE / name).exists()]
    if missing:
        raise SystemExit(f"Missing prepared input(s) in {SOURCE}: {', '.join(missing)}")

    for name in PREPARED_INPUTS:
        shutil.copy2(SOURCE / name, OUT / name)

    land_geometry = base.load_geojson(OUT / "sodermalm_osm_island_epsg3006.geojson")["features"][0]["geometry"]
    shoreline = base.polygon_boundary_rings(land_geometry)
    all_shore_points = [p for ring in shoreline for p in ring]
    px = np.array([p[0] for p in all_shore_points])
    py = np.array([p[1] for p in all_shore_points])
    center = (float((px.min() + px.max()) * 0.5), float((py.min() + py.max()) * 0.5))
    radius = float(max(math.hypot(x - center[0], y - center[1]) for x, y in all_shore_points) + base.BUFFER_M)

    samples, water_level = base.terrain_samples(OUT / "sodermalm_cylinder_contours.geojson", center)
    geo = OUT / f"sodermalm_cylinder_{TAG}_quad_disk.geo"
    msh = OUT / f"sodermalm_cylinder_{TAG}_quad_disk.msh"
    mesh_path = OUT / f"sodermalm_cylinder_{TAG}.nmsh"

    base.write_gmsh_geo(geo, radius, shoreline, center)
    base.run([base.GMSH, "-2", str(geo), "-format", "msh2", "-o", str(msh)])
    nodes, quads, lines = base.parse_msh(msh)
    mesh_info = base.write_nmsh(mesh_path, nodes, quads, lines, samples, water_level, shoreline, center, radius)

    building_stl_info = base.write_sodermalm_building_stl(
        OUT / f"sodermalm_buildings_{TAG}.stl",
        OUT / "sodermalm_buildings.geojson",
        samples,
        water_level,
        shoreline,
        center,
        radius,
    )

    if not skip_figures:
        base.render_mesh(FIG / f"sodermalm_cylinder_{TAG}_mesh_topdown.png", nodes, quads, radius, shoreline, center)
        base.render_mesh(FIG / f"sodermalm_cylinder_{TAG}_mesh_topdown_notext.png", nodes, quads, radius, shoreline, center)
        base.render_terrain_mesh_3d(
            FIG / f"sodermalm_cylinder_{TAG}_mesh_3d.png",
            nodes,
            quads,
            lines,
            samples,
            water_level,
            shoreline,
            center,
            radius,
        )
        base.render_terrain_only(
            FIG / f"sodermalm_terrain_height_contours_{TAG}_notext.png",
            shoreline,
            OUT / "sodermalm_cylinder_contours.geojson",
            center,
            radius,
            samples,
            water_level,
            annotate=False,
        )
        base.render_3d_scene(
            FIG / f"sodermalm_cylinder_{TAG}_3d_buildings.png",
            shoreline,
            OUT / "sodermalm_buildings.geojson",
            center,
            radius,
            samples,
            water_level,
        )
    checker_output = base.run_mesh_checker(mesh_path)

    metadata = {
        "target": "2x finer horizontal local p2 run, same vertical setup",
        "center_epsg3006": center,
        "radius_m": radius,
        "buffer_m": base.BUFFER_M,
        "domain_height_above_lowest_m": base.DOMAIN_HEIGHT_ABOVE_LOWEST_M,
        "vertical_stretch": base.VERTICAL_STRETCH,
        "nz": base.NZ,
        "land_mesh_size_m": base.LAND_MESH_SIZE_M,
        "effective_land_p2_xy_spacing_m": base.LAND_MESH_SIZE_M / 2.0,
        "water_mesh_size_m": base.WATER_MESH_SIZE_M,
        "effective_water_p2_xy_spacing_m": base.WATER_MESH_SIZE_M / 2.0,
        "shoreline_simplify_m": base.SHORELINE_SIMPLIFY_M,
        "inflow_from_degrees": base.INFLOW_FROM_DEG,
        "inflow_arc_width_degrees": base.INFLOW_ARC_WIDTH_DEG,
        "outflow_arc_width_degrees": base.OUTFLOW_ARC_WIDTH_DEG,
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
        "mesh_checker_summary": checker_output,
    }
    (OUT / f"sodermalm_cylinder_{TAG}_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
