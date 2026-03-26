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

**File transfer node:** use `dardel-ftn` for rsync. Use `dardel` for sbatch/monitoring.
Account: `naiss2025-3-39`, partition: `main`, 128 cores/node.

**Cluster scripts** (in `cluster/`):
- `build_neko_channel.sh` — SLURM build job (sbatch wrapper)
- `sync_to_dardel.sh` — rsync source → Dardel
- `sync_from_dardel.sh <run_name> [--all]` — rsync results → egidius
- `job_template.sh` — Dardel run job template (copy and edit CASE → actual name)

## Running the two-phase channel example

```bash
cd examples/two_phase_channel
# Generate mesh (once)
genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 81 18 27 .true. .false. .true.
mv box.nmsh box_phys_81x18x27.nmsh
source ../../setup-env-channel.sh --egidius   # or --local
makeneko turb_channel_two_phase.f90
mpirun -np 4 ./neko turb_channel_two_phase_v4.case
```

For production runs on egidius, run from a simulations directory:
```bash
mkdir -p /lscratch/sieburgh/simulations/channel_test_v4
cd /lscratch/sieburgh/simulations/channel_test_v4
cp ~/code/neko-multiphase-channel/examples/two_phase_channel/*.case .
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
makeneko turb_channel_single_phase.f90
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
| `examples/two_phase_channel/turb_channel_two_phase.f90` | **Phase 1** user module (original, no n̂ fix). Use for Phase 1 reference cases. |
| `examples/two_phase_channel/turb_channel_two_phase_p2.f90` | **Phase 2** user module: identical to Phase 1 except adds a second GS averaging pass on n̂ before `div()`. The diff is documented in the file header. Compile with `makeneko turb_channel_two_phase_p2.f90` for all `p2_*` cases. |
| `examples/two_phase_channel/turb_channel_two_phase_laminar.case` | Laminar + We=1 (σ=0.3): blew up at t=0.90 TU (same CSF instability as turbulent cases) |
| `examples/two_phase_channel/turb_channel_two_phase_we1.case` | Phase 1 turbulent + We=1 (σ=0.3): blew up; superseded by p2_we1 |
| `examples/two_phase_channel/turb_channel_two_phase_we10.case` | Phase 1 turbulent + We=10 (σ=0.03): blew up at t=20.44 (Δt/Δt_cap=0.69); superseded by p2_we10 |
| `examples/two_phase_channel/turb_channel_two_phase_sigma0.case` | Phase 1 σ=0 CDI-only quality test; κ_rms→64; root cause confirmed as element-face C0 kink in n̂ |
| `examples/two_phase_channel/turb_channel_two_phase_p2_sigma0.case` | **Phase 2** σ=0 fix-verification: 108×18×36 mesh, ε=0.09, R=0.4. Target: κ_rms≈5.0 stable |
| `examples/two_phase_channel/turb_channel_two_phase_p2_we10.case` | **Phase 2** We=10 (σ=0.04): first σ>0 test with Fortran fix; 108×18×36 mesh |
| `examples/two_phase_channel/turb_channel_two_phase_p2_we1.case` | **Phase 2** We=1 (σ=0.4): primary production case; 108×18×36 mesh |
| `examples/two_phase_channel/turb_channel_two_phase_restart.case` | We=1.33 restart from `fluid00004.chkp` (t=20→25), R=0.4, y_c=0 (centre); Phase 1 blow-up reference |
| `examples/two_phase_channel/turb_channel_two_phase_restart_off.case` | We=1.33 restart, R=0.4, y_c=0.3 (off-centre, log-law region) |
| `examples/two_phase_channel/turb_channel_two_phase_v4.case` | v4: We=730 high-We reference (completed, Phase 1) |
| `examples/two_phase_channel/turb_channel_single_phase.case` | Single-phase spin-up to t=25, checkpoints every 5 TU. **Must be re-run on 108×18×36 mesh for Phase 2.** |
| `examples/two_phase_channel/turb_channel_single_phase.f90` | Fluid-only user module for single-phase spin-up |
| `examples/two_phase_channel/postprocess_single_phase.py` | Single-phase postprocessing: ekin plot + mean velocity profile |
| `examples/two_phase_channel/animate_blowup.py` | Animation: φ/κ/\|u\| panels. Flags: `--stride N`, `--mesh p1\|p2`, `--kappa-scale`. |
| `examples/turb_channel/turb_channel.f90` | Reference: channel IC source |
| `examples/spurious_currents_multiphase/spurious_currents.f90` | Reference: CSF/CDI source |
| `cluster/build_neko_channel.sh` | Dardel SLURM build job |
| `cluster/sync_to_dardel.sh` | rsync source → Dardel |
| `cluster/sync_from_dardel.sh` | rsync results ← Dardel |
| `cluster/job_template.sh` | Dardel run job template |
| `setup-env-channel.sh` | Environment setup (local/egidius/cluster) |
| `build-neko-channel.sh` | Build script |

## Phase 2 workflow (new mesh + Fortran fix)

### Mesh

```bash
# 108x18x36 — uniform xz (Δx=Δz=0.1164), 3.1 elements across interface
genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 108 18 36 .true. .false. .true.
mv box.nmsh box_phys_108x18x36.nmsh
```

Key parameters for Phase 2 cases: `ε=0.09`, `R=0.4`, `γ=0.05`.
- 4ε/Δ = 3.1 elements across interface (vs 1.8 in Phase 1)
- We=10: σ=0.04, Δt_cap≈0.079, Δt/Δt_cap≈0.051
- We=1: σ=0.4, Δt_cap≈0.025, Δt/Δt_cap≈0.16

### Phase 2 strategy

**Phase 2 investigates the finer mesh only — the Fortran fix is deferred.**
The original `turb_channel_two_phase.f90` is used for all Phase 2 cases.
`turb_channel_two_phase_p2.f90` (with extra GS pass on n̂) exists for later
comparison but is not used yet. The question is: does better interface
coverage (3.1 elements vs 1.8) alone change the blow-up behaviour?

### Compile for Phase 2

```bash
cd examples/two_phase_channel
source ../../setup-env-channel.sh --egidius   # or --cluster on Dardel
makeneko turb_channel_two_phase.f90           # original code, no Fortran fix
```

### Verification sequence

1. Run `p2_sigma0.case` on Dardel (1 TU, 128 ranks). Observe κ_rms — expect still ~64 if element-face artifact dominates mesh resolution, or lower if coverage helps.
2. If verified, sync to Dardel, build, re-run single-phase spin-up on 108×18×36 mesh (128 ranks),
   then run p2_we10 and p2_we1.

### Single-phase spin-up on new mesh (Dardel)

```bash
# Update turb_channel_single_phase.case: change mesh_file to box_phys_108x18x36.nmsh
mkdir $SCRATCH_DIR/channel_p2_single_phase
cd    $SCRATCH_DIR/channel_p2_single_phase
cp $SRC/turb_channel_single_phase.case .
cp $SRC/box_phys_108x18x36.nmsh .
cp $SRC/neko .
srun -u -n 128 ./neko turb_channel_single_phase.case
# fluid00004.chkp written at t=20; copy to each Phase 2 run directory
```

**MPI rank count for restart must match spin-up.** New-mesh spin-up uses 128 ranks →
all Phase 2 restart cases use `srun -n 128`.

## Style and architecture

Same as the parent Neko project — see the Fortran style rules in the
upstream CONTRIBUTING.md. Key reminder: Cray compiler on Dardel enforces
`implicit none` strictly; all pointer variables must be explicitly declared
in every subroutine.
