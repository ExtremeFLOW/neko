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

Current production mesh geometry:

```sh
SODERMALM_OUTPUT_SUFFIX=p7_xy35_inlet220_nz10 \
SODERMALM_TARGET_LABEL="production mesh geometry" \
SODERMALM_LAND_MESH_SIZE_M=35 \
SODERMALM_WATER_MESH_SIZE_M=225 \
SODERMALM_INFLOW_ARC_WIDTH_DEG=220 \
SODERMALM_NZ=10 \
SODERMALM_DOMAIN_HEIGHT_M=1000 \
SODERMALM_SKIP_FIGURES=1 \
python3 build_sodermalm_cylinder_2x_xy.py
```

This gives `35 m` land elements and a `1 km` domain height. The qualified
diagnostic uses polynomial order `3`, or about `11.7 m` effective x-y spacing
over land. The same geometry gives about `5 m` effective spacing at order `7`.

## Build Brinkman Cache On LUMI

The LUMI run uses a cached building indicator field to avoid repeating the STL
distance search in every solver job.

```sh
python3 make_template_field_from_nmsh.py \
  generated_p7_xy35_inlet220_nz10/sodermalm_cylinder_p7_xy35_inlet220_nz10.nmsh \
  cache_p3_xy35_inlet220_nz10_allbuildings_centered_sw50_h65/template.f00000 \
  --polynomial-order 3

python3 make_sodermalm_building_cache.py \
  --generated generated_p7_xy35_inlet220_nz10 \
  --template-field cache_p3_xy35_inlet220_nz10_allbuildings_centered_sw50_h65/template.f00000 \
  --cache-prefix cache_p3_xy35_inlet220_nz10_allbuildings_centered_sw50_h65/bldp3_ctr50 \
  --smooth-width 50 \
  --transition-placement centered \
  --max-height 65
```

By default the cache builder keeps every land-valid building footprint. Do not
apply area/edge filters unless a specific mesh-resolution diagnostic proves
that unresolved footprint fragments are causing instability. `run_neko.sh`
builds this cache automatically when it is absent.

Verify that every flat building base intersects the terrain before running:

```sh
python3 check_building_ground_contact.py generated_p7_xy35_inlet220_nz10
python3 render_building_contact_closeups.py generated_p7_xy35_inlet220_nz10
```

The audit must report zero buildings with an air gap. The close-ups show the
most and least embedded bases against the terrain surface.

## Run Neko Locally

```sh
cd run_sodermalm_wind_2x_xy_t10_np6_soft_buildings
./run_neko.sh
```

Optional overrides:

```sh
NP=4 MPIEXEC=mpirun NEKO_BIN=/path/to/neko ./run_neko.sh
```

The case runs to `t = 100`, uses the Högdalen vertical wind profile as both
initial and inlet velocity, applies it on the untapered 220-degree inlet, and
writes fields every 10 seconds under `fields`.

## Run Baseline On LUMI

```sh
sbatch run_lumi_sodermalm_baseline.sbatch
```

The baseline Slurm script expects the case under
the LUMI scratch case directory. Do not commit the site-local Slurm script,
because it contains the project account.

Current qualified baseline (LUMI job `21977386`):

- mesh: `generated_p7_xy35_inlet220_nz10/sodermalm_cylinder_p7_xy35_inlet220_nz10.nmsh`
- building object: cached centered all-land-valid building indicator
- polynomial order: `3`
- dealiasing: enabled
- wind: full Högdalen/SW profile as initial condition and inlet condition
- inlet: `220 deg`, no taper, centered on SW wind-from direction
- side arcs: `outflow+dong`
- top: `normal_outflow`
- bottom: no-slip terrain/water bottom
- Reynolds number: `100`
- Brinkman penalty: `10`
- Brinkman limits: `[0, 65]`
- Brinkman ramp time: `30 s`
- Brinkman PDE filter radius: `20 m`
- wind ramp: disabled
- `target_cfl`: `0.15`
- `max_timestep`: `0.02`
- `max_dt_increase_factor`: `1.1`
- output interval: `10 s`
- checkpoint interval: `25 s`

This clean `t = 0` to `100 s` run completed past the former repeatable
`t ~= 57 s` failure. Its final maximum velocity was `13.37 m/s`, divergence
RMS was `9.35e-3`, and pressure-residual RMS was `8.92e-3`.

Do not use `normal_outflow` on the cylindrical side arcs with the default
non-full-stress formulation. It selects Neko's axis-aligned tangential
constraint, which is not valid on a curved boundary. Building-generated
perturbations exposed this as a repeatable backflow instability near
`t = 57 s`. Use the energy-stable `outflow+dong` condition on both open side
arcs; `normal_outflow` remains appropriate on the horizontal top.

LUMI cached-Brinkman runs need enough host memory and the current device-MPI
build. Keep the site-local Slurm wrapper and project account outside Git.

Do not enable `SODERMALM_RAMP_TIME` with this top-boundary setup unless the
top/open reference velocity is ramped consistently too. Ramping only the
initial/inlet wind introduced a mismatch at the top/open boundary in testing.

## Render The Heatmap

```sh
GENERATED_PATH=generated_p7_xy35_inlet220_nz10 \
FIELD_PATH=run_sodermalm_wind_2x_xy_t10_np6_soft_buildings/fields/field0.f00010 \
GEOMETRY_FIELD=run_sodermalm_wind_2x_xy_t10_np6_soft_buildings/fields/field0.f00000 \
LIFT_M=10 SAMPLE_ONLY=1 \
NPZ_OUT=run_sodermalm_wind_2x_xy_t10_np6_soft_buildings/renders/t100_terrain_plus10.npz \
python3 sample_terrain_relative_velocity_npz.py

VMIN=0 BUILDING_ALPHA=95 \
python3 render_velocity_field_heatmap.py \
  generated_p7_xy35_inlet220_nz10 \
  run_sodermalm_wind_2x_xy_t10_np6_soft_buildings/renders/t100_terrain_plus10.npz \
  run_sodermalm_wind_2x_xy_t10_np6_soft_buildings/renders/t100_terrain_plus10.png
```

The first script samples the solution at `terrain + 10 m`. The second renders
those samples as a continuous field, then overlays the shoreline, buildings,
inlet arc, and color scale. Use one shared uncapped maximum across frames when
making a movie so colors remain comparable in time.
