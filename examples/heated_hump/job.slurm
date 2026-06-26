#!/usr/bin/env bash
#SBATCH --job-name=heated_hump
#SBATCH --nodes=8                 # EDIT: e.g. 8 nodes
#SBATCH --ntasks-per-node=64      # EDIT: cores per node on your cluster
#SBATCH --time=24:00:00
#SBATCH --partition=cpu           # EDIT: your CPU partition
#SBATCH --account=YOUR_ACCOUNT    # EDIT
#SBATCH --output=hump_%j.log
#
# Flow over a wall-mounted hump, low-Mach (CPU-only), recycling turbulent inflow,
# heated wall (x>11), gravity. ~512 MPI ranks for the 53k-elem / order-7 mesh.

set -euo pipefail
JSONFORTRAN="${JSONFORTRAN:-$HOME/json-fortran-install}"
NEKO_DIR="${NEKO_DIR:-$HOME/neko-lowmach}"

# --- load the SAME toolchain you built with (EDIT) ---
# module load gcc openmpi openblas
export LD_LIBRARY_PATH="$JSONFORTRAN/lib:${LD_LIBRARY_PATH:-}"

cd "$SLURM_SUBMIT_DIR"

# build the driver if not already built on this filesystem
[ -x ./neko ] || "$NEKO_DIR/install/bin/makeneko" heated_hump.f90

# optional one-time sanity check (cheap):
# srun -n 8 ./neko smoke.case

srun ./neko heated_hump.case      # srun uses $SLURM_NTASKS ranks automatically
