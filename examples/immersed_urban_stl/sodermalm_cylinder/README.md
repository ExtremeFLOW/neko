# Sodermalm Cylinder: P9 Sharp Mask

The retained configuration covers Sodermalm, Gamla Stan, and Riddarholmen.
Terrain forms the mesh bottom; buildings form the Brinkman solid. This is a
production-development case, not a claim of validated building-scale accuracy.
The retained baseline completed a p7-to-p9 continuation from t = 225 to 500 s.
This is not evidence that a fresh p9 start or different geometry has been tested.

## Current Settings

- Conforming quadrilateral disk extruded into hexes; radius about 2.54 km.
- Land elements 35 m, water elements 225 m; p9 gives nominal 3.89 m land spacing.
  GLL spacing is nonuniform.
- Original mesh: 347,160 hexes (34,716 base quads).
- Ten vertically stretched layers; top 350 m above the lowest terrain/water.
- Wind from SW (225 degrees); full height-dependent profile at initialization
  and inlet: `4.467 * (max(z, 1) / 51)^0.16` m/s. No wind ramp or inlet taper.
- Inlet: 180 degrees, bearings 135 through 315. Remaining cylindrical side
  zones: `outflow+dong`; top: `normal_outflow`; bottom: no-slip.
- All land-valid building footprints, no area/edge filtering. Heights are
  limited to 2-65 m; bases extend 2 m below sampled footprint minimum terrain.
- Sharp raw mask (zero distance-transition width), PDE filter radius 10 m.
  Explicit Brinkman: penalty 10, limits [0, 65], forcing ramp 30 s.
- HPFRT: two modes, weight 2; Re 100; dealiasing enabled.
- CFL target 0.15; initial dt 0.001, maximum dt 0.02, growth factor 1.1.
- Restart to t = 500 s; fields every 10 s, checkpoints every 25 s.

The definitive solver parameters are in `run/sharp_mask.case`.
The successful run used a mesh named `p7_xy35_inlet220_nz10`, but inspection of
its binary boundary labels confirms **180 degrees**, not 220. The defaults
here preserve that actual configuration. The generator still splits circle
curves at exact zone endpoints and checks the binary labels; the 220-degree
option remains supported, but is not the completed p9 baseline.

## Prepare

Requirements: Python 3.11+ with NumPy and Pillow; Matplotlib for rendering.
For new geometry: Gmsh, `mesh_checker`,
and `rsvg-convert` on PATH (or set `GMSH`, `MESH_CHECKER`, `RSVG`).
Neko must include the current Brinkman filter/ramp and HPFRT support.
FFmpeg (on PATH or via `imageio-ffmpeg`) is needed only for movies.

### Reproduce the Successful Continuation

Reuse the original mesh and its prepared GeoJSON/metadata, not a regenerated
mesh with the same dimensions. Place or symlink these in `generated_sharp_mask/`,
with the mesh named `sodermalm_cylinder_sharp_mask.nmsh`. Its SHA-256 is:

```text
846b42f0fb4c8de0d196388e7f2dce77a7cd9b1850b5a7a25f5e8a7893c0098d
```

Symlink the original p7 checkpoint at t = 225.0001350978 to
`run/restart.chkp`. The tested checkpoint has 347,160 elements, lx = 8, and
22,751,477,944 bytes. Neko interpolates the checkpoint's fields and histories
to p9 on the same element mesh; do not specify `restart_mesh_file`.
The building mask must be rebuilt at p9 nodes, not interpolated from p7:

```sh
python3 prepare_case.py
```

Preparation verifies the actual inlet, original mesh hash, checkpoint header
and completeness, then creates the p9 cache and records its input hashes.
`--check` also verifies the cache order and size without regenerating it.
These structural checks do not replace a numerical audit of a new checkpoint.
Private run locations and allocation scripts remain outside the repository.

### Generate New Geometry

This is for a new run, not for restarting the successful checkpoint. Remove
`restart_file` from the case to start from the full wind profile at t = 0.

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
Preparation creates the p9 building cache in `cache_sharp_mask/` and records
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
LIFT_M=20 NPZ_OUT=renders/terrain_plus20.npz \
python3 sample_terrain_relative_velocity_npz.py

python3 render_velocity_field_heatmap.py generated_sharp_mask \
  renders/terrain_plus20.npz renders/terrain_plus20.png
```

For all selected outputs, with a single uncapped scale across frames:

```sh
python3 render_movie.py run/fields/field0.f????? \
  --geometry-field run/fields/field0.f00000 \
  --lift 20 --fps 5 --out renders/terrain_plus20_movie
```

Use `--lift 30` for terrain + 30 m. The movie directory must be new.
Outputs include PNG frames, sampled NPZ files, and `velocity.mp4`.
Rendering retains the latest 32-band Cool-to-Warm colormap with logarithmic
mapping (reference speed 0.025 m/s), black building footprints with open
courtyards, and a labeled colorbar. No titles or streamlines. Building interiors
are hidden only in the visualization, not modified in the solver outputs.
Five frames/s matches the latest fast movie; use `--fps 1` for slower playback.
Add `--gif` for the same frames as a GIF; `--width` controls video/GIF size
(default 1400 px), without reducing the saved PNGs.
The scale is shared across all frames and never caps the sampled/displayed
maximum. An explicitly supplied `VMAX` below that maximum is rejected.
The inlet marker uses the verified 180 degrees; for a different mesh, set
`INLET_WIDTH_DEG` to its verified width rather than trusting old metadata.

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
