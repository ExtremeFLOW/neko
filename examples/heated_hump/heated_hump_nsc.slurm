#!/bin/bash
#SBATCH -N 16                       # 53k elems @ order 7 = 27.5M pts; 16 nodes(512 ranks)~54k pts/rank. 8-32 all fine.
#SBATCH -t 1000
#SBATCH -A liu-compute-2026-4
#SBATCH --exclusive
#SBATCH -J heated_hump
#SBATCH -o heated_hump_%j.log
#SBATCH -e heated_hump_error_%j.log
echo "=========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_JOB_NUM_NODES   Start: $(date)"
echo "=========================================="

module load buildenv-gcc/2022a-eb 2>/dev/null

# json-fortran runtime (same as your bfs setup)
export LD_LIBRARY_PATH=$HOME/json-fortran-install/jsonfortran-gnu-9.2.0/lib:$LD_LIBRARY_PATH

# IMPORTANT: use the LOW-MACH fork install, NOT $HOME/neko-install (that one has no lowmach scheme)
export PATH=$HOME/neko-lowmach-install/bin:$PATH

cd /proj/liu2/users/x_adigh/heated_hump/     # <-- EDIT: dir with heated_hump.{f90,case}, hump.nmsh
echo "Working directory: $(pwd)"

NPROCS=$(( SLURM_JOB_NUM_NODES * 32 ))        # 32 cores/node on Tetralith
echo "Running with $NPROCS MPI ranks"

# build the case-specific executable if not already built (makeneko -> ./neko)
if [ ! -x ./neko ]; then
    echo "compiling driver with makeneko ..."
    makeneko heated_hump.f90 || { echo "makeneko failed"; exit 1; }
fi

echo "Starting Neko (./neko, the makeneko-compiled binary with the user module)..."
mpirun -np $NPROCS ./neko heated_hump.case    # NB: ./neko  -- NOT the generic 'neko'

rc=$?
echo "=========================================="
[ $rc -eq 0 ] && echo "Neko completed OK. End: $(date)" || { echo "Neko FAILED. End: $(date)"; exit 1; }
echo "=========================================="
