# Immersed boundary buildings in a Gmsh terrain box

This example is a smoke test for a GIS-to-Neko urban immersed-boundary workflow.

- The outer fluid mesh is a coarse terrain-following box generated as a Gmsh
  HEX27 mesh and converted through `gmsh2nek` and `rea2nbin`.
- The buildings are generated as a separate STL and used as a Brinkman
  `boundary_mesh` immersed boundary.

The goal is only to verify STL ingestion, Brinkman mask construction, and a
short coarse run. This is not a production urban-flow setup.

## Local Geometry Package

The minimal geometry material needed to regenerate the example lives under:

```sh
geometry_pipeline/
```

Included inputs:

- `geometry_pipeline/input/steep_pilot_layers.gpkg`: clipped pilot GIS package.

Included scripts:

- `geometry_pipeline/scripts/build_urban_stl.py`: creates `terrain.stl`,
  `buildings.stl`, and related STL diagnostics from the local GeoPackage.
- `geometry_pipeline/scripts/make_neko_gmsh_mesh.py`: creates the coarse Gmsh
  HEX27 terrain-following mesh.
- `geometry_pipeline/generate_geometry.sh`: wrapper that runs the STL and mesh
  generation steps and links the results into this example folder.

Generated geometry is written to `geometry_pipeline/generated/`, which is
ignored by git.

## Generate Mesh And STL

From this directory:

```sh
./geometry_pipeline/generate_geometry.sh
```

This creates a coarse `16 x 16 x 5` terrain-following mesh. Edit the constants
near the top of the Python scripts to change the geometry or mesh resolution.

The wrapper expects the usual GIS/mesh tools on `PATH`: `ogr2ogr`, `gdal_grid`,
`gmsh2nek`, `rea2nbin`, and `mesh_checker`. If the Neko tools are installed in
`/usr/local/neko_install/bin`, add that directory to `PATH` before running.

The wrapper creates or updates these local symlinks:

- `urban_terrain.nmsh`
- `urban_terrain_surface.stl`
- `urban_buildings.stl`

To relink/check already generated artifacts without rebuilding them:

```sh
./prepare.sh
```

## Run Neko

After the geometry is prepared:

```sh
neko immersed_urban_stl.case
```

Successful output should show that the building STL was read as a
`boundary_mesh`, the Brinkman source was initialized, and fields were written to
`fields/`.

## Render The Reference Image

After a Neko run has produced `fields/field0.f00000`, render the terrain, solid
depth-buffered buildings, and continuous near-ground streamlines:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python3 render_urban_depth_scene.py
```

The renderer is intentionally self-contained and uses a simple z-buffer so the
building STL is drawn as solid 3D geometry instead of Matplotlib's approximate
3D polygon sorting.
