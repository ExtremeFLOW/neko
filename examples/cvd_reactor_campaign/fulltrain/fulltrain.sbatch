#!/bin/bash
#SBATCH -N 10
#SBATCH -t 1000
#SBATCH -A liu-compute-2026-4
#SBATCH --exclusive
#SBATCH -J fulltrain
#SBATCH -o fulltrain_%j.log
#SBATCH -e fulltrain_error_%j.log
# 10 nodes x 32 = 320 ranks for 68,376 elements ~ 214 elements/rank -- good
# parallel efficiency for this element count (23.4M GLL at N=6, ~73k pts/rank).
#
# At N=6 dt is CFL-limited to ~7e-6 by the 0.63 mm side-wall BL elements, so
# one 16 h job advances t by ~0.17 = 0.19 of a flow-through (0.908). Reaching
# 3 flow-throughs is a ~16-job campaign (~34,000 core-h). Each follow-on job
# needs "restart_file": "fluid000NN.chkp" added to the case.
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
# Add the LOW-MACH Neko fork to PATH (NOT $HOME/neko-install -- no lowmach).
# $HOME/neko-lowmach/install is the Jul 21 build, commit aa3daf8b: the same
# code heated_hump runs, and the one whose log prints the Scalar block before
# Fluid/Pressure/Momentum. The stale Jun 26 build was parked at
# $HOME/neko-lowmach-install.stale-jun26 and the old name is now a symlink
# here, so either path resolves to this build.
export PATH=$HOME/neko-lowmach/install/bin:$PATH
export LD_LIBRARY_PATH=$HOME/neko-lowmach/install/lib:$LD_LIBRARY_PATH
# Change to project directory
cd /proj/liu-compute-2026-4/users/x_adigh/neko_low_mach_solver/fulltrain/
echo "Working directory: $(pwd)"
# Get total number of MPI ranks
NPROCS=$(( SLURM_JOB_NUM_NODES * 32 ))
echo "Running with $NPROCS MPI ranks on $SLURM_JOB_NUM_NODES nodes"
echo "=========================================="
echo "Starting Neko simulation..."
echo "=========================================="
# Run Neko -- use ./neko (the makeneko-compiled binary with the user module,
# i.e. the buoyancy source term and the parabolic velocity inlet),
# NOT the generic 'neko'.
mpirun -np $NPROCS ./neko fulltrain.case
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
