#!/bin/sh
#SBATCH -A naiss2025-3-39
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -J advecting_drop_test
#SBATCH -p main

### MODULES & ENV ####
module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p
######################

# Set up Neko environment
export NEKO_DIR=/cfs/klemming/projects/supr/kthmech/eriksie/software/neko-cpu
export PATH=$NEKO_DIR/bin:$PATH
export LD_LIBRARY_PATH=$NEKO_DIR/lib:$LD_LIBRARY_PATH

### RUN ###
mkdir -p "visualization_output"

makeneko advecting_drop.f90
chmod +x neko
srun -u -n 64 ./neko advecting_drop.case