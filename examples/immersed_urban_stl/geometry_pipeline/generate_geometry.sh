#!/usr/bin/env bash

set -e

HERE=$(cd "$(dirname "$0")" && pwd)
EXAMPLE=$(cd "$HERE/.." && pwd)
PYTHON=/opt/homebrew/Caskroom/miniconda/base/bin/python3
PATH=/usr/local/neko_install/bin:/usr/local/bin:/opt/homebrew/bin:$PATH

"$PYTHON" "$HERE/scripts/build_urban_stl.py"
"$PYTHON" "$HERE/scripts/make_neko_gmsh_mesh.py"

cd "$HERE/generated"
mesh_checker urban_terrain.nmsh

cd "$EXAMPLE"
mkdir -p fields
rm -f urban_terrain.nmsh urban_terrain_surface.stl urban_buildings.stl
ln -sf geometry_pipeline/generated/urban_terrain.nmsh urban_terrain.nmsh
ln -sf geometry_pipeline/generated/terrain.stl urban_terrain_surface.stl
ln -sf geometry_pipeline/generated/buildings.stl urban_buildings.stl
