# Södermalm Cylinder Wind Smoke Case

This folder contains the local Södermalm immersed-boundary wind example used for
quick CPU tests and GPU production-pilot runs. The cutout includes Södermalm
plus the upstream Gamla Stan island group. The mesh is a cylindrical extruded
hex mesh with a terrain-following bottom, open top, and building geometry
provided as an STL or cached Brinkman indicator field.

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

Default local mesh target:

- land element size: `110 m`, giving about `55 m` effective p2 x-y spacing
- water element size: `225 m`
- vertical elements: `5`
- domain top: `350 m` above the lowest terrain/water level
- polynomial order in the case: `2`
- default wind direction: SW, using meteorological convention
  (`225 deg`, north is `0/360 deg`, clockwise positive)
- inlet arc: `180 deg`, centered on the wind direction

Current p7 production-pilot mesh target:

```sh
SODERMALM_OUTPUT_SUFFIX=p7_xy35_inlet220 \
SODERMALM_TARGET_LABEL="p7 production-pilot mesh" \
SODERMALM_LAND_MESH_SIZE_M=35 \
SODERMALM_WATER_MESH_SIZE_M=225 \
SODERMALM_INFLOW_ARC_WIDTH_DEG=220 \
SODERMALM_NZ=5 \
SODERMALM_DOMAIN_HEIGHT_M=350 \
SODERMALM_SKIP_FIGURES=1 \
python3 build_sodermalm_cylinder_2x_xy.py
```

This gives about `5 m` effective x-y spacing over land at polynomial order `7`.

## Build Brinkman Cache On LUMI

The LUMI run uses a cached building indicator field to avoid repeating the STL
distance search in every solver job.

```sh
python3 make_sodermalm_building_cache.py \
  --generated generated_2x_xy \
  --template-field fields/field0.f00000 \
  --cache-prefix fields/sodermalm_buildings_2x_allbuildings_sw50_h65_distance_cache \
  --smooth-width 50
```

By default the cache builder keeps every land-valid building footprint. Do not
apply area/edge filters unless a specific mesh-resolution diagnostic proves
that unresolved footprint fragments are causing instability.

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
the LUMI scratch case directory. Do not commit the site-local Slurm script,
because it contains the project account.

Current stable p7 production-pilot setup:

- mesh: `generated_p7_xy35_inlet220/sodermalm_cylinder_p7_xy35_inlet220.nmsh`
- building object: cached all-land-valid building Brinkman indicator
- polynomial order: `7`
- dealiasing: enabled
- wind: full Högdalen/SW profile as initial condition and inlet condition
- inlet: `220 deg`, no taper, centered on SW wind-from direction
- side arcs: `outflow+dong`
- top: `normal_outflow`
- bottom: no-slip terrain/water bottom
- Reynolds number: `100`
- Brinkman penalty: `0.1`
- Brinkman limits: `[0, 65]`
- Brinkman ramp time: `300 s`
- wind ramp: disabled
- `target_cfl`: `0.05`
- `max_timestep`: `0.01`
- `max_dt_increase_factor`: `1.05`
- output interval: `10 s`
- checkpoint interval: `50 s`

Do not enable `SODERMALM_RAMP_TIME` with this top-boundary setup unless the
top/open reference velocity is ramped consistently too. Ramping only the
initial/inlet wind introduced a mismatch at the top/open boundary in testing.

Useful override pattern:

```sh
RUN_NAME=lumi_p7_allbuild_t500_cfl005_br300_nowindramp \
MESH=../generated_p7_xy35_inlet220/sodermalm_cylinder_p7_xy35_inlet220.nmsh \
CACHE_NAME=bldp7_all0 \
BRINKMAN_OBJECT=cache \
END_TIME=500 OUTPUT_VALUE=10 OUTPUT_CHECKPOINTS=true \
CHECKPOINT_CONTROL=simulationtime CHECKPOINT_VALUE=50 \
POLYNOMIAL_ORDER=7 DEALIAS=true \
SODERMALM_INITIAL_WIND_SCALE=1 SODERMALM_RAMP_TIME=0 \
SODERMALM_INLET_ARC_WIDTH_DEG=220 SODERMALM_INLET_TAPER_DEG=0 \
BRINKMAN_LIMIT_MAX=65 BRINKMAN_PENALTY=0.1 BRINKMAN_RAMP_TIME=300 \
TARGET_CFL=0.05 MAX_TIMESTEP=0.01 MAX_DT_INCREASE_FACTOR=1.05 \
  sbatch run_lumi_sodermalm_baseline.sbatch
```

The p7 `t=500` pilot job with these settings passed the previous instability
window at `t ~= 60` without solver errors.

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
