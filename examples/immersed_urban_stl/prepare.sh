#!/usr/bin/env bash

set -e

HERE=$(cd "$(dirname "$0")" && pwd)
GENERATED="$HERE/geometry_pipeline/generated"

cd "$HERE"
mkdir -p fields

rm -f urban_terrain.nmsh urban_terrain_surface.stl urban_buildings.stl
ln -sf geometry_pipeline/generated/urban_terrain.nmsh urban_terrain.nmsh
ln -sf geometry_pipeline/generated/terrain.stl urban_terrain_surface.stl
ln -sf geometry_pipeline/generated/buildings.stl urban_buildings.stl

cd "$GENERATED"
mesh_checker urban_terrain.nmsh
