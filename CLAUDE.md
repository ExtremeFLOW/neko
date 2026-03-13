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

See the Dardel workflow in the original `spurious-currents` CLAUDE.md.
Use `--cluster` flag. Neko installs to:
`/cfs/klemming/projects/supr/kthmech/eriksie/builds/neko-channel/`

## Running the two-phase channel example

```bash
cd examples/two_phase_channel
cp ../turb_channel/box.nmsh .      # mesh is gitignored
source ../../setup-env-channel.sh --egidius   # or --local
makeneko turb_channel_two_phase.f90
mpirun -np 4 ./neko turb_channel_two_phase.case
```

For production runs on egidius, run from a simulations directory:
```bash
mkdir -p /lscratch/sieburgh/simulations/channel_test
cd /lscratch/sieburgh/simulations/channel_test
cp ~/code/neko-multiphase-channel/examples/two_phase_channel/*.case .
cp ~/code/neko-multiphase-channel/examples/two_phase_channel/neko .
cp ~/code/neko-multiphase-channel/examples/turb_channel/box.nmsh .
mpirun -np 32 ./neko turb_channel_two_phase.case
```

## Key source files

| Path | Purpose |
|------|---------|
| `examples/two_phase_channel/turb_channel_two_phase.f90` | User module (main file) |
| `examples/two_phase_channel/turb_channel_two_phase.case` | Local test case (N=5, end_time=5) |
| `examples/turb_channel/turb_channel.f90` | Reference: channel IC source |
| `examples/spurious_currents_multiphase/spurious_currents.f90` | Reference: CSF/CDI source |
| `setup-env-channel.sh` | Environment setup (local/egidius/cluster) |
| `build-neko-channel.sh` | Build script |

## Style and architecture

Same as the parent Neko project — see the Fortran style rules in the
upstream CONTRIBUTING.md. Key reminder: Cray compiler on Dardel enforces
`implicit none` strictly; all pointer variables must be explicitly declared
in every subroutine.
