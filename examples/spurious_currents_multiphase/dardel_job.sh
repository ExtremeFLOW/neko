#!/bin/bash
# Slurm job template for spurious currents sigma sweep on Dardel.
# This file is instantiated by run_sigma_sweep.sh --cluster,
# which replaces the SIGMA and END_TIME placeholders before submitting.
#
# Edit the #SBATCH directives below to match your project and time needs
# before running the sweep.

#SBATCH -A <YOUR_ACCOUNT>       # e.g. kthmech
#SBATCH -p gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4     # 4 MI250X GPUs per node on Dardel
#SBATCH -t 04:00:00             # adjust per sigma (sigma=0.01 needs ~24h)
#SBATCH -J sc_SIGMA
#SBATCH -o sc_SIGMA.log
#SBATCH -e sc_SIGMA.log

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

export NEKO_SPURIOUS_PREFIX=$KTHMECH_PROJECT/software/neko-spurious

echo "Starting sigma=SIGMA  end_time=END_TIME  at $(date)"
srun "$NEKO_SPURIOUS_PREFIX/bin/neko" spurious_currents.case
echo "Done at $(date)"
