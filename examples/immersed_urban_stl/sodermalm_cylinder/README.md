# Södermalm Cylinder Wind Smoke Case

This folder contains the local Södermalm immersed-boundary wind example used for
quick CPU tests. The cutout includes Södermalm plus the upstream Gamla Stan
island group. The mesh is a cylindrical extruded hex mesh with a 180-degree
inlet arc, terrain-following bottom, open top, and building geometry
provided as an STL for the Brinkman source term.

## Prepare GIS Cutout

```sh
python3 expand_sodermalm_with_gamla_stan.py
```

The script expands the prepared Södermalm land mask with Gamla Stan and
Riddarholmen, then reclips buildings and terrain contours from the local
GeoPackage. It keeps a
`sodermalm_only_epsg3006.geojson` base mask so repeated runs do not accumulate
duplicate island geometry. A small close/open is applied to merge the very tight
old-town gaps that are below this mesh resolution.

## Generate Geometry

```sh
python3 build_sodermalm_cylinder_2x_xy.py
```

The script reads prepared GIS cutouts from `input/prepared`, writes the mesh and
building STL to `generated_2x_xy`, produces diagnostic figures in
`figures_2x_xy`, and runs `mesh_checker` on the generated `.nmsh`.

Current mesh target:

- land element size: `110 m`, giving about `55 m` effective p2 x-y spacing
- water element size: `225 m`
- vertical elements: `5`
- domain top: `350 m` above the lowest terrain/water level
- polynomial order in the case: `2`
- default wind direction: SW, using meteorological convention
  (`225 deg`, north is `0/360 deg`, clockwise positive)
- inlet arc: `180 deg`, centered on the wind direction

## Build Brinkman Cache On LUMI

The LUMI run uses a cached building indicator field to avoid repeating the STL
distance search in every solver job.

```sh
python3 make_sodermalm_building_cache.py \
  --generated generated_2x_xy \
  --field fields/field0.f00000 \
  --out fields/sodermalm_buildings_2x_clean50_area25_distance_cache0.f00000 \
  --smooth-step 50 \
  --min-footprint-area 25 \
  --min-edge 2 \
  --simplify 2
```

## Run Neko Locally

```sh
cd run_sodermalm_wind_2x_xy_t10_np6_soft_buildings
./run_neko.sh
```

Optional overrides:

```sh
NP=4 MPIEXEC=mpirun NEKO_BIN=/path/to/neko ./run_neko.sh
```

The case currently runs to `t = 100`,
uses the Högdalen vertical wind profile, applies it on the 180-degree inlet,
uses a soft Brinkman building mask, and writes fields under `fields`.

## Run Baseline On LUMI

```sh
sbatch run_lumi_sodermalm_baseline.sbatch
```

The baseline Slurm script expects the case under
`/scratch/project_465002721/daotuana/neko-sodermalm-hip/case`, uses the
Neko HIP build at `install/neko-hip-1998-device-mpi-o0`, reuses the cached
`generated_2x_xy` Brinkman mask, and runs the successful diagnostic setup:

- one LUMI GPU node, four MPI ranks
- prefilled Högdalen/SW wind profile initial condition
- end time `t = 25`
- output every `5` time units
- Brinkman limits `[0, 50]`, penalty `0.1`

Useful overrides:

```sh
END_TIME=100 OUTPUT_VALUE=10 RUN_NAME=run_lumi_sodermalm_t100 \
  sbatch run_lumi_sodermalm_baseline.sbatch
```

## Render The Heatmap

```sh
cd run_sodermalm_wind_2x_xy_t10_np6_soft_buildings
python3 render_terrain_plus10_heatmap.py
```

The renderer samples the final velocity field at `terrain + 10 m`, overlays the
Södermalm shoreline, buildings, and inlet arc, and writes
`renders/sodermalm_velocity_terrain_plus10.png`.

The diagnostic plot that exposed the terrain/building imprint used a focused
velocity range:

```sh
FIELD_PATH=fields/field0.f00006 \
GEOMETRY_FIELD=fields/field0.f00000 \
SHOW_STREAMLINES=0 \
VMIN=2.7 VMAX=5.1 \
BUILDING_ALPHA=50 \
python3 render_terrain_plus10_heatmap.py
```
