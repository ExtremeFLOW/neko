# Immersed Urban STL

## Requirements

Put the local GIS input here:

```sh
geometry_pipeline/input/steep_pilot_layers.gpkg
```

Required commands on `PATH`:

```sh
ogr2ogr
gdal_grid
gmsh2nek
rea2nbin
mesh_checker
neko
```

## Generate Geometry

From this directory:

```sh
./geometry_pipeline/generate_geometry.sh
./prepare.sh
```

This creates:

```sh
urban_terrain.nmsh
urban_terrain_surface.stl
urban_buildings.stl
```

## Generate Brinkman Mask

First create the template field:

```sh
cd fine_3x/run_mask
rm -rf fields
mkdir -p fields
neko immersed_urban_stl.case
cd ../..
```

Then create the footprint-based mask:

```sh
python3 geometry_pipeline/scripts/make_brinkman_mask.py
```

The mask is written to:

```sh
fine_3x/upload/urban_brinkman_mask0.f00000
fine_3x/upload/urban_brinkman_mask0.nek5000
```

## Run Neko

For the default local case:

```sh
neko immersed_urban_stl.case
```

For the `fine_3x` mesh with the generated mask, copy the files from
`fine_3x/upload/` into the run directory and use the Brinkman object:

```json
{
  "type": "file",
  "file_name": "fields/urban_brinkman_mask0.fld",
  "field_name": "s01"
}
```

## Render

After `fields/field0.f00000` exists:

```sh
python3 render_urban_depth_scene.py
```

Output:

```sh
urban_depth_full_ground_streamlines_view_restored.png
```
