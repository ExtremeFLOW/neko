#!/bin/bash
# Slurm job template for spurious currents sigma sweep on Dardel.
# This file is instantiated by run_sigma_sweep.sh --cluster,
# which replaces the SIGMA and END_TIME placeholders before submitting.

#SBATCH -A naiss2025-3-39
#SBATCH -p main
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128   # CPU cores per node on Dardel
#SBATCH -t 04:00:00             # adjust per sigma (sigma=0.01 needs ~24h)
#SBATCH -J sc_SIGMA
#SBATCH -o sc_SIGMA.log
#SBATCH -e sc_SIGMA.log

PROJECT_DIR=/cfs/klemming/projects/supr/kthmech/eriksie
NEKO_SPURIOUS_PREFIX=$PROJECT_DIR/builds/neko-spurious

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

export PATH="$NEKO_SPURIOUS_PREFIX/bin:$PATH"

echo "Compiling spurious_currents.f90 at $(date)"
makeneko spurious_currents.f90

echo "Starting sigma=SIGMA  end_time=END_TIME  at $(date)"
srun ./neko spurious_currents.case
echo "Done at $(date)"
