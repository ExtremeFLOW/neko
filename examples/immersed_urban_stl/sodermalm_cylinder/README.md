# Sodermalm Cylinder: Sharp Mask

The retained configuration covers Sodermalm, Gamla Stan, and Riddarholmen.
Terrain forms the mesh bottom; buildings form the Brinkman solid. This is a
production-development case, not a claim of validated building-scale accuracy.

## Current Settings

- Conforming quadrilateral disk extruded into hexes; radius about 2.54 km.
- Land elements 35 m, water elements 225 m; p7 gives nominal 5 m land spacing.
  GLL spacing is nonuniform.
- Ten vertically stretched layers; top 350 m above the lowest terrain/water.
- Wind from SW (225 degrees); full height-dependent profile at initialization
  and inlet: `4.467 * (max(z, 1) / 51)^0.16` m/s. No wind ramp or inlet taper.
- Inlet: 220 degrees, bearings 115 through 335. Remaining cylindrical side
  zones: `outflow+dong`; top: `normal_outflow`; bottom: no-slip.
- All land-valid building footprints, no area/edge filtering. Heights are
  limited to 2-65 m; bases extend 2 m below sampled footprint minimum terrain.
- Sharp raw mask (zero distance-transition width), PDE filter radius 10 m.
  Explicit Brinkman: penalty 10, limits [0, 65], forcing ramp 30 s.
- HPFRT: two modes, weight 2; Re 100; dealiasing enabled.
- CFL target 0.15; initial dt 0.001, maximum dt 0.02, growth factor 1.1.
- Fresh start to t = 500 s; fields every 10 s, checkpoints every 25 s.

The definitive solver parameters are in `run/sharp_mask.case`.
The inlet fix is geometric: circle curves are split at exact zone endpoints,
and the binary mesh is checked after writing. Older files named `inlet220`
can actually have a 180-degree inlet. Do not reuse their cache or checkpoint
with regenerated geometry. The current long sharp-mask run predates this
geometry correction; this directory combines its solver settings with the fix.

## Prepare

Requirements: Python 3.11+ with NumPy and Pillow; Gmsh, `mesh_checker`,
and `rsvg-convert` on PATH (or set `GMSH`, `MESH_CHECKER`, `RSVG`).
Neko must include the current Brinkman filter/ramp and HPFRT support.
FFmpeg is needed only for movies.

Preserve these prepared GIS inputs in `input/prepared/` (EPSG:3006):

- `sodermalm_osm_island_epsg3006.geojson`
- `sodermalm_osm_water_cut_epsg3006.geojson`
- `sodermalm_buildings.geojson`
- `sodermalm_cylinder_contours.geojson`

From this directory:

```sh
python3 generate_geometry.py
python3 check_building_ground_contact.py generated_sharp_mask
python3 prepare_case.py
```

Generation writes `generated_sharp_mask/` and geometry figures to
`figures_sharp_mask/`. The contact audit must find no air gaps.
Preparation creates the p7 building cache in `cache_sharp_mask/` and records
its input hashes. It rejects unverified existing caches. To omit geometry
figures, use `SODERMALM_SKIP_FIGURES=1`.

## Build and Run

Build on the appropriate host before submitting a solver job:

```sh
cd run
makeneko sodermalm_wind.f90
bash run.sh mpirun -np 6
```

The launch above illustrates a local MPI command, not a production allocation.
On a cluster, pass the existing known-working launcher/binding arguments to
`run.sh` from a private scheduler script. The script performs input checks,
does not compile, and refuses to overwrite previous fields. Keep account IDs
and scheduler wrappers outside Git.

The user callback logs maximum velocity/location and divergence diagnostics
every simulated second. Set `SODERMALM_DIAGNOSTIC_INTERVAL` to change this
interval. No solver or boundary experiment is enabled by an environment flag.

## Heatmaps and Movies

Render existing outputs on the cluster where they reside; download only media.
For one frame:

```sh
FIELD_PATH=run/fields/field0.f00001 \
GEOMETRY_FIELD=run/fields/field0.f00000 \
LIFT_M=10 NPZ_OUT=renders/terrain_plus10.npz \
python3 sample_terrain_relative_velocity_npz.py

python3 render_velocity_field_heatmap.py generated_sharp_mask \
  renders/terrain_plus10.npz renders/terrain_plus10.png
```

For all selected outputs, with a single uncapped scale across frames:

```sh
python3 render_movie.py run/fields/field0.f????? \
  --geometry-field run/fields/field0.f00000 \
  --lift 10 --out renders/terrain_plus10_movie
```

Use `--lift 30` for terrain + 30 m. The movie directory must be new.
Outputs include PNG frames, sampled NPZ files, and `velocity.mp4`.
Rendering retains the logarithmic urban-flow colormap, building/shoreline
overlays and colorbar, without titles or streamlines. Building interiors are
not masked out.

Sampling is approximate: it selects nearest-height GLL values in horizontal
bins (about 12 m), applies the existing display smoothing, and renders a
continuous heatmap. It is not spectral-element interpolation. The color range
covers the sampled/displayed slice, not the whole 3D field; use solver
diagnostics for global maxima.

## Checks

The small regression tests protect the corrected boundary zoning and retained
case settings; they do not launch simulations.

```sh
python3 -m unittest discover -s tests -v
python3 prepare_case.py --check
```

Production notes, job logs, GIS data, caches, results, and scheduler scripts
stay outside version control.
