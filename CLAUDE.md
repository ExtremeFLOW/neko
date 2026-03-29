# CLAUDE.md — neko-multiphase-channel

This is branch `eriksie/multiphase/two-phase-channel` of the Neko fork at
`git@github.com:ExtremeFLOW/neko.git`. It builds on the validated CSF
implementation from `eriksie/multiphase/spurious-currents` (curvature sign fix,
commit 637fc27) and adds a 3D turbulent two-phase channel flow case.

## Build system

Neko uses Autotools. Channel-specific scripts handle environment + build:

```bash
# Set up environment (pick one)
source setup-env-channel.sh --local     # macOS with gfortran/clang
source setup-env-channel.sh --egidius   # KTH Mechanics workstation
source setup-env-channel.sh --cluster   # Dardel HPC (Cray + ROCm)

# Build
./build-neko-channel.sh --local
./build-neko-channel.sh --egidius
./build-neko-channel.sh --egidius --clean   # full rebuild
```

## Machines

### Local (macOS)

Install prefix: `$HOME/local/neko-channel`
json-fortran: `$HOME/local` (built separately)

### Egidius (KTH Mechanics workstation)

**Storage:** NFS home (`~`) is slow. Everything lives on local XFS scratch:

| Purpose | Path |
|---------|------|
| Source code | `/lscratch/sieburgh/code/neko-multiphase-channel/` |
| Built Neko | `/lscratch/sieburgh/local/neko-channel/` |
| json-fortran | `/lscratch/sieburgh/local/jsonfortran-gnu-9.2.1/` |
| Simulation runs | `/lscratch/sieburgh/simulations/<run_name>/` |

`~/code`, `~/local`, `~/data`, `~/simulations` are all symlinks to `/lscratch/sieburgh/`.

**Compiler:** gfortran 14.2.0, OpenMPI 5.0.7, 32 cores.
No job scheduler — run `mpirun -np <N>` directly.

**json-fortran** was built from source with cmake:
```bash
cd ~/local/src/json-fortran && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/lscratch/sieburgh/local -DSKIP_DOC_GEN=ON
make -j$(nproc) install
# installs to ~/local/jsonfortran-gnu-9.2.1/
```

**Known issue:** cmake leaves `libdir` and `includedir` empty in the generated
`json-fortran.pc`. The `.pc` file must be patched manually after install:
```
libdir=/lscratch/sieburgh/local/jsonfortran-gnu-9.2.1/lib
includedir=/lscratch/sieburgh/local/jsonfortran-gnu-9.2.1/lib
```
(already done; verify with `pkg-config --cflags json-fortran` after sourcing the env)

### Dardel (production HPC)

Use `--cluster` flag. Key paths:

| Purpose | Path |
|---------|------|
| Source code | `$KTHMECH_PROJECT/src/neko-multiphase-channel/` |
| Built Neko | `$KTHMECH_PROJECT/builds/neko-channel/` |
| Simulation runs | `$SCRATCH_DIR/<run_name>/` |

Where `KTHMECH_PROJECT=/cfs/klemming/projects/supr/kthmech/eriksie` and
`SCRATCH_DIR=/cfs/klemming/scratch/e/eriksie`.

**One-time setup:**
```bash
cd $KTHMECH_PROJECT/src
git clone git@github.com:ExtremeFLOW/neko.git neko-multiphase-channel
cd neko-multiphase-channel && git checkout eriksie/multiphase/two-phase-channel
sbatch $KTHMECH_PROJECT/scripts/build_neko_channel.sh
```

**Node roles — use the right node for each task:**

| Task | Node | Command |
|------|------|---------|
| File transfers (rsync, scp) | `dardel-ftn` | `rsync ... dardel-ftn:...` |
| Job submission, monitoring | `dardel` | `sbatch`, `squeue`, `git pull` |
| Builds, compilation | Compute node via SLURM | `sbatch build_neko_channel.sh` |
| Simulations | Compute node via SLURM | `sbatch job_*.sh` |
| Mesh generation | **egidius** (generate + rsync) | `genmeshbox ... && rsync ... dardel-ftn:...` |

**Do not run genmeshbox or makeneko interactively on the Dardel login node.**
Generate meshes on egidius and transfer via `dardel-ftn`.

Account: `naiss2025-3-39`, partition: `main`, 128 cores/node.

**Cluster scripts** (in `cluster/`):
- `build_neko_channel.sh` — SLURM build job for Neko framework (makeneko + libs)
- `build_neko_two_phase.sh` — SLURM build job for two-phase user module → `neko_two_phase`
- `sync_to_dardel.sh` — rsync source → Dardel
- `sync_from_dardel.sh <run_name> [--all]` — rsync results → egidius
- `job_template.sh` — Dardel run job template (copy and edit CASE → actual name)
- `job_channel_p2_single_phase.sh` — L1 (108×18×36) single-phase spin-up
- `job_channel_p2_sigma0.sh` — L1 σ=0, ε=0.09 (informative; not convergence series)
- `job_channel_p2_sigma0_eps053.sh` — L1 σ=0, ε=0.053 (convergence series point)
- `job_channel_p2_we10.sh` — L1 We=10 two-phase
- `job_channel_p2_we1.sh` — L1 We=1 two-phase
- `job_channel_p3_single_phase.sh` — L2 (144×24×48) single-phase spin-up
- `job_channel_p3_sigma0.sh` — L2 σ=0 CDI diagnostic, ε=0.04 (full 5 TU)
- `job_channel_p3_sigma0_p2.sh` — L2 σ=0 _p2 variant (NO-OP for σ=0; not useful)
- `job_channel_p3_we10.sh` — L2 We=10 two-phase
- `job_channel_p3_we1.sh` — L2 We=1 two-phase
- `job_channel_l3_single_phase.sh` — L3 (192×32×64) single-phase spin-up, 4 nodes/512 ranks
- `job_channel_l3_sigma0.sh` — L3 σ=0 CDI diagnostic, ε=0.03 (full 5 TU)
- `job_channel_l3_we10.sh` — L3 We=10 two-phase
- `job_channel_l3_we1.sh` — L3 We=1 two-phase

## Running the two-phase channel example

```bash
cd examples/two_phase_channel
# Generate mesh (once)
genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 81 18 27 .true. .false. .true.
mv box.nmsh box_phys_81x18x27.nmsh
source ../../setup-env-channel.sh --egidius   # or --local
makeneko src/turb_channel_two_phase.f90
mpirun -np 4 ./neko turb_channel_two_phase_v4.case
```

For production runs on egidius, run from a simulations directory:
```bash
mkdir -p /lscratch/sieburgh/simulations/channel_test_v4
cd /lscratch/sieburgh/simulations/channel_test_v4
cp ~/code/neko-multiphase-channel/examples/two_phase_channel/cases/81x18x27/*.case .
cp ~/code/neko-multiphase-channel/examples/two_phase_channel/neko .
cp ~/code/neko-multiphase-channel/examples/two_phase_channel/box_phys_81x18x27.nmsh .
mpirun -np 32 ./neko turb_channel_two_phase_v4.case
```

### Turbulent restart (two-step) — standard workflow for all turbulent two-phase cases

**All turbulent two-phase cases restart from `fluid00004.chkp`.** This gives genuinely
developed turbulence at t=20 rather than a synthetic Reichardt IC. The laminar case
is the only exception (it uses a Poiseuille IC with no perturbations).

```bash
# Step 1: single-phase spin-up — generates fluid00004.chkp at t=20
makeneko src/turb_channel_single_phase.f90
mpirun -np 16 ./neko turb_channel_single_phase.case
# Verify: u_max in ekin.csv fluctuating ~1.35–1.45 with no trend at t=20
# (mean profile matches Reichardt Re_tau=180 — confirmed via postprocess_single_phase.py)

# Step 2: two-phase restart — reads fluid00004.chkp, injects drop analytically
# Replace <CASE> with restart, we10, we1, restart_off, etc.
mkdir /lscratch/sieburgh/simulations/channel_test_<CASE>
cd /lscratch/sieburgh/simulations/channel_test_<CASE>
cp ~/code/neko-multiphase-channel/examples/two_phase_channel/{turb_channel_two_phase_<CASE>.case,neko,box_phys_81x18x27.nmsh} .
cp /lscratch/sieburgh/simulations/channel_single_phase/fluid00004.chkp .
mpirun -np 16 ./neko turb_channel_two_phase_<CASE>.case
```

**Restart mechanics:** Neko restores t=20 from the checkpoint. `end_time: 25.0` is
absolute — the simulation runs from t=20 to t=25. The fluid IC is skipped (velocity
from checkpoint); the scalar IC places the drop analytically. Verify in neko.log:
`Restarting from checkpoint ... Time: 20.0` and `T : [20.0, 25.0]`.

**MPI rank count for restart must match spin-up.** The single-phase run used 16 ranks;
fluid00004.chkp is partitioned for 16 ranks. Restart with `mpirun -np 16`.

**Checkpoint compatibility fix (in this fork):** `src/io/chkp_file.f90` and
`src/case.f90` were patched to allow restarting from a fluid-only checkpoint into a
scalar case. If `have_scalar=0` in the checkpoint, the scalar IC is applied from user
code instead of erroring. The `chkp_t` struct carries a `scalar_was_read` flag.

**drop_center_y parameter:** `case.scalar.drop_center_y` (default 0.0) sets the
wall-normal offset of the drop centre. Used by the restart_off case (y_c=0.3).

## Key source files

| Path | Purpose |
|------|---------|
| `examples/two_phase_channel/src/turb_channel_two_phase.f90` | Original user module — used for all convergence series runs |
| `examples/two_phase_channel/src/turb_channel_two_phase_p2.f90` | Variant with extra GS on n̂ — **NO-OP**: face nodes already consistent after GS(∇φ)+normalize; kink is intra-element, unreachable by GS |
| `examples/two_phase_channel/src/turb_channel_single_phase.f90` | Fluid-only user module for single-phase spin-up (all meshes) |
| `examples/two_phase_channel/cases/81x18x27/` | All case files for the baseline mesh (single-phase + two-phase runs) |
| `examples/two_phase_channel/cases/108x18x36/` | Case files for L1 mesh (108×18×36, ε=0.053) |
| `examples/two_phase_channel/cases/144x24x48/` | Case files for L2 mesh (144×24×48, ε=0.04) |
| `examples/two_phase_channel/cases/192x32x64/` | Case files for L3 mesh (192×32×64, ε=0.03) |
| `examples/two_phase_channel/postprocess/postprocess_single_phase.py` | Single-phase postprocessing: ekin plot + mean velocity profile |
| `examples/two_phase_channel/postprocess/animate_blowup.py` | Animation: φ/κ/\|u\| panels. Flags: `--stride N`, `--mesh p1\|p2\|p3`, `--kappa-scale` |
| `examples/two_phase_channel/postprocess/postprocess_sigma0.py` | σ=0 diagnostics: κ_rms/φ time-series + field snapshots. Args: `--run`, `--R`, `--eps`, `--mesh`, `--no-snapshots` |
| `examples/two_phase_channel/figures/` | Output figures and animations (gitignored, generated locally) |
| `examples/turb_channel/turb_channel.f90` | Reference: channel IC source |
| `examples/spurious_currents_multiphase/spurious_currents.f90` | Reference: CSF/CDI source |
| `examples/turb_channel/turb_channel.f90` | Reference: channel IC source |
| `examples/spurious_currents_multiphase/spurious_currents.f90` | Reference: CSF/CDI source |
| `cluster/build_neko_channel.sh` | Dardel SLURM build job |
| `cluster/build_neko_two_phase.sh` | Dardel build: `turb_channel_two_phase.f90` → `neko_two_phase` |
| `cluster/build_neko_two_phase_p2.sh` | Dardel build: `turb_channel_two_phase_p2.f90` → `neko_two_phase_p2` |
| `cluster/sync_to_dardel.sh` | rsync source → Dardel |
| `cluster/sync_from_dardel.sh` | rsync results ← Dardel |
| `cluster/job_template.sh` | Dardel run job template |
| `setup-env-channel.sh` | Environment setup (local/egidius/cluster) |
| `build-neko-channel.sh` | Build script |

## Convergence series — L1/L2/L3/L4 workflow

Isotropic meshes with constant ε/Δxz=0.457, ε/Δ_GLL=3.2 (N=7 fixed). R=0.4, γ=0.05.
All cases use `turb_channel_two_phase.f90`.

| Level | Mesh | ε | Δxz | R/ε | Elements | Nodes (est.) |
|-------|------|---|-----|-----|----------|--------------|
| L1 | 108×18×36 | 0.053 | 0.1164 | 7.5 | 70k | 1 |
| L2 | 144×24×48 | 0.040 | 0.0873 | 10 | 166k | 2 |
| L3 | 192×32×64 | 0.030 | 0.0654 | 13 | 393k | 4 |
| L4 | 288×48×96 | 0.020 | 0.0436 | 20 | 1327k | ~14 |

**L4 note:** The 4/3 refinement factor breaks at L4 (64×4/3=85.3, not integer). L4 uses
factor 3/2 from L3 (nz: 64→96). ε=0.02 is a clean round number. L4 is planned for
high-resolution visualisation and validation, not yet submitted. Will require ~14 nodes
(1792 ranks); may need to go higher still for publication-quality visualisations.

### L1 (108×18×36) — spin-up COMPLETED; sigma0 cases COMPLETED (2026-03-27/28)

Key parameters: ε=0.053 (convergence series) or ε=0.09 (informative only), R=0.4, γ=0.05.
Spin-up: job 18985538 completed. `fluid00004.chkp` in `$SCRATCH_DIR/channel_p2_single_phase/`.

Two-phase jobs (jobs 19003682–19003685) submitted 2026-03-28. Results synced to egidius:
- `channel_p2_sigma0_eps053`: COMPLETED t=20→25, κ_rms≈56 (ε=0.053, convergence series point)
- `channel_p2_sigma0`: COMPLETED t=20→25, κ_rms≈50 (ε=0.09, informative)
- `channel_p2_we10`: synced
- `channel_p2_we1`: synced

**MPI rank count for restart must match spin-up.** L1 spin-up: 128 ranks → all L1 restarts use `srun -n 128`.

### Build two-phase binary on Dardel

```bash
bash cluster/sync_to_dardel.sh
ssh dardel "cp $KTHMECH_PROJECT/src/neko-multiphase-channel/cluster/build_neko_two_phase.sh \
              $KTHMECH_PROJECT/scripts/ && \
            sbatch $KTHMECH_PROJECT/scripts/build_neko_two_phase.sh"
```

Note: `build_neko_two_phase_p2.sh` also exists but the _p2 variant is a no-op (see `turb_channel_two_phase_p2.f90` header). Not needed for current runs.

### L2 (144×24×48) — spin-up COMPLETED (job 19002586, 2026-03-28); TURBULENT validated

`box_phys_144x24x48.nmsh` generated on egidius 2026-03-27. Spin-up: 2 nodes / 256 ranks.
u_max (last 5) = 1.37 ± 0.01 — TURBULENT. `fluid00004.chkp` in `$SCRATCH_DIR/channel_p3_single_phase/`.
Results synced to `/lscratch/sieburgh/simulations/channel_p3_single_phase/`.

**Submit L2 two-phase (sigma0 pair + we10 + we1):**
```bash
bash cluster/sync_to_dardel.sh
ssh dardel "cp $KTHMECH_PROJECT/src/neko-multiphase-channel/cluster/job_channel_p3_*.sh \
              $KTHMECH_PROJECT/scripts/"
ssh dardel "sbatch $KTHMECH_PROJECT/scripts/job_channel_p3_sigma0.sh"
ssh dardel "sbatch $KTHMECH_PROJECT/scripts/job_channel_p3_sigma0_p2.sh"
ssh dardel "sbatch $KTHMECH_PROJECT/scripts/job_channel_p3_we10.sh"
ssh dardel "sbatch $KTHMECH_PROJECT/scripts/job_channel_p3_we1.sh"
```

**MPI rank count:** L2 spin-up uses 256 ranks (2 nodes) → all L2 restarts use `srun -n 256`.
Note: 1 node/128 ranks OOM at 1296 elem/rank; 2 nodes/256 ranks = 648 elem/rank.

### L3 (192×32×64) — spin-up COMPLETED (job 19003203, 2026-03-28); TURBULENT validated

`box_phys_192x32x64.nmsh` generated on egidius 2026-03-27. Spin-up: 4 nodes / 512 ranks.
u_max (last 5) = 1.34 ± 0.00 — TURBULENT. `fluid00004.chkp` in `$SCRATCH_DIR/channel_l3_single_phase/`.
Results synced to `/lscratch/sieburgh/simulations/channel_l3_single_phase/`.

**Submit L3 two-phase (sigma0 + we10 + we1):**
```bash
bash cluster/sync_to_dardel.sh
ssh dardel "cp $KTHMECH_PROJECT/src/neko-multiphase-channel/cluster/job_channel_l3_*.sh \
              $KTHMECH_PROJECT/scripts/"
ssh dardel "sbatch $KTHMECH_PROJECT/scripts/job_channel_l3_sigma0.sh"
ssh dardel "sbatch $KTHMECH_PROJECT/scripts/job_channel_l3_we10.sh"
ssh dardel "sbatch $KTHMECH_PROJECT/scripts/job_channel_l3_we1.sh"
```

**MPI rank count:** L3 spin-up uses 512 ranks (4 nodes) → all L3 restarts use `srun -n 512`.
Note: 2 nodes/256 ranks OOM at 1536 elem/rank; 4 nodes/512 ranks = 768 elem/rank.

## Style and architecture

Same as the parent Neko project — see the Fortran style rules in the
upstream CONTRIBUTING.md. Key reminder: Cray compiler on Dardel enforces
`implicit none` strictly; all pointer variables must be explicitly declared
in every subroutine.
