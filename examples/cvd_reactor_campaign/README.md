# SiC CVD reactor campaign — meshes and setups

Retrieval bundle for the low-Mach SiC CVD simulations (H2 carrier, 200 mbar,
1800 K walls), August 2026. Gas-phase chemistry lives in a separate private
repo (`sic_cvd_chem`); this directory holds the flow side: meshes, case
setups, mesh generation, and the restart-chain tooling.

## meshes/cvd_meshes.zip

| mesh | elements | for |
|---|---|---|
| `hotzone_trip_sem.nmsh` | 29,616 hex27 | hot zone + inlet extension + floor trip hump (zones: 1 inlet, 2 outlet, 3 wallCold, 4 wallHot, 5 substrate) |
| `train.nmsh` | 68,376 hex27 | full train: inlet pipe + constricted hump + O-grid diffusor + hot zone + extended outlet (minSJ 0.732) |

## hotzone_trip/ — hot-zone cases (Re_H 2000, T_hot 6 = 1800 K)

`hotzone_trip.{f90,case}` — flow + heated walls (spatial + temporal
smoothstep ramps; see `SETUP_NOTES.md` for the ramp-sizing rules and the
measured failures each ramp prevents).
`hotzone_species.{f90,case}` — + 5 passive species (STEP 1 of the chemistry
track; reacting versions in the private chem repo).

## fulltrain/ — full-train turbulent-jet case (Re_pipe 4000)

- `fulltrain.f90` + `fulltrain_run1.case` — fresh start, N = 6, low-Mach H2,
  both wall ramps, gravity on w.
- `fulltrain_run2_restart.case` — the tuned chained-restart configuration:
  `restart_file`, projection_space_size 20, `max_timestep 2.2e-5` (dt pinned
  under the cap = fixed-dt behavior with a CFL safety valve),
  checkpoint_value 0.05.
- `fulltrain.sbatch` (10 nodes) / `fulltrain2.sbatch` (14 nodes, 24 h).
- `handover.sh` — login-node watcher: waits for a checkpoint to finish
  writing (size-stable check), cancels the old job, submits the next.
- `slice_train.py` / `slice_train_z.py` — element-aware y=0 side view and
  z=0 top view slicers (pysemtools).

KNOWN BUG for restarts: with `variable_timestep: true`, a checkpoint restart
skips all simulationtime outputs in the first 0.1 time units (placeholder
dt = 1.0 enters set_counter's 0.1*dt guard). Keep checkpoint_value well
below 0.1. Full analysis in the neko_bug notes.

## mesh_tools/ — full-train mesh generation chain

`ogrid_cvd.py` (geometry + O-grid butterfly morphing, hump constriction,
hot-zone cross-section tables) -> `ogrid_smooth_cvd.py` (PSI smoothing,
morph-blended) -> `build_train.py` (151 stations, SMID diffusor
distribution, zone retag, substrate staircase) + `tune_diffusor.py`
(streamwise-shear quality tuning) -> `build_nmsh.py` (hex27 curving BEFORE
mm->m scaling -> .nmsh). `preview_train.py` renders the geometry;
`build_diffusor.py` / `build_full_train.py` are earlier stages kept for
reference.

## Cluster layout (sigma, account liu-compute-2026-4)

`/proj/liu-compute-2026-4/users/x_adigh/neko_low_mach_solver/{fulltrain,
full_train_run_2, sic_skeletal, ...}`; build: makeneko from
`$HOME/neko-lowmach/install` (this branch, Jul-21 build), see the sbatch
files for the module/paths pattern.
