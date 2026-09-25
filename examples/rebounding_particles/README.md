# Rebounding inertial particles

This example demonstrates Lagrangian particle tracking with inertial particles,
nonlinear drag, and elastic rebounds from a wall boundary. The fluid is frozen
and initially at rest, so the particle motion can be compared against the
one-dimensional analytical drag solution used by the plotting scripts.

The folder contains two equivalent case files:

- `rebounding_particles_csv.case` reads the initial particles from
  `particles.csv` and writes particle trajectories as CSV.
- `rebounding_particles_h5.case` defines the initial particles directly in the
  case file and writes particle trajectories as HDF5.

## Running the cases

To run the CSV input/output path, first create the input particle file and then
launch Neko:

```bash
python3 particles_init.py
neko rebounding_particles_csv.case
python3 trajectory_plot_csv.py
```

To run the JSON input and HDF5 output path:

```bash
neko rebounding_particles_h5.case
python3 trajectory_plot_h5.py
```

The postprocessing scripts generate `particle_trajectories.png` and
`velocity.png`.

## LPT output options

The LPT output file is configured with a base filename and a format:

```json
"output_filename": "tracers",
"output_format": "h5"
```

If `output_format` is omitted, it defaults to `csv`. Supported formats are
`csv`, `h5`, and `hdf5`.

The optional `snapshots_per_file` setting controls how many written snapshots
are stored in each output file:

```json
"snapshots_per_file": "all"
```

keeps the current schema with all snapshots in one output sequence, while

```json
"snapshots_per_file": 1
```

writes each output snapshot to a separate file. Larger positive integers group
that many output snapshots per file. The `output_control` and `output_value`
settings still control when snapshots are written.

## Mesh

The provided `box.nmsh` mesh is a small box and wall
zones suitable for the elastic rebound demonstration.
