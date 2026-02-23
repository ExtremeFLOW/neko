#!/bin/bash
# Slurm job for oscillating droplet on Dardel (CPU partition).
# Compiles the user file and runs the simulation.
# Output is written to the working directory (should be on scratch).

#SBATCH -A naiss2025-3-39
#SBATCH -p main
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH -t 04:00:00
#SBATCH -J osc_drop
#SBATCH -o osc_drop.log
#SBATCH -e osc_drop.log

PROJECT_DIR=/cfs/klemming/projects/supr/kthmech/eriksie
NEKO_SPURIOUS_PREFIX=$PROJECT_DIR/builds/neko-spurious

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

export PATH="$NEKO_SPURIOUS_PREFIX/bin:$PATH"

echo "Compiling oscillating_droplet.f90 at $(date)"
makeneko oscillating_droplet.f90

echo "Starting oscillating droplet simulation at $(date)"
srun ./neko oscillating_droplet.case
echo "Done at $(date)"
