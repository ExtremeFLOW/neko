#!/bin/bash
#SBATCH -N 45
#SBATCH -t 1000
#SBATCH -A liu-compute-2026-4
#SBATCH --exclusive
#SBATCH -J heated_hump
#SBATCH -o heated_hump_%j.log
#SBATCH -e heated_hump_error_%j.log
echo "=========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Running on node(s): $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "Number of tasks: $SLURM_NTASKS"
echo "Start time: $(date)"
echo "=========================================="
# Load the build environment (needed for MPI runtime)
module load buildenv-gcc/2022a-eb 2>/dev/null
# Set library paths
export LD_LIBRARY_PATH=$HOME/json-fortran-install/jsonfortran-gnu-9.2.0/lib:$LD_LIBRARY_PATH
# Add the LOW-MACH Neko fork to PATH  (NOT $HOME/neko-install -- that one has no lowmach)
export PATH=$HOME/neko-lowmach-install/bin:$PATH
# Change to project directory
cd /proj/liu-compute-2026-4/users/x_adigh/neko_low_mach_solver/
echo "Working directory: $(pwd)"
# Get total number of MPI ranks
NPROCS=$(( SLURM_JOB_NUM_NODES * 32 ))
echo "Running with $NPROCS MPI ranks on $SLURM_JOB_NUM_NODES nodes"
echo "=========================================="
echo "Starting Neko simulation..."
echo "=========================================="
# Run Neko -- use ./neko (the makeneko-compiled binary with the user module),
# NOT the generic 'neko' (which has no recycling/heating/gravity hooks).
mpirun -np $NPROCS ./neko heated_hump.case
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "Neko simulation completed successfully!"
    echo "End time: $(date)"
    echo "=========================================="
else
    echo "=========================================="
    echo "Neko simulation failed!"
    echo "End time: $(date)"
    echo "=========================================="
    exit 1
fi
