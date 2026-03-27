#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -J build_neko_two_phase
#SBATCH -p main
#SBATCH -o /cfs/klemming/projects/supr/kthmech/eriksie/logs/%J_build_neko_two_phase.log

# Compile the two-phase user module on Dardel.
# Run AFTER build_neko_channel.sh has completed (makeneko must be installed).
#
# Output: $SRC/neko_two_phase  (used by all p2/p3 two-phase job scripts)
#
# Submit from the Dardel login node:
#   sbatch $KTHMECH_PROJECT/scripts/build_neko_two_phase.sh

set -e

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

SRC=$KTHMECH_PROJECT/src/neko-multiphase-channel/examples/two_phase_channel

echo "Compiling turb_channel_two_phase.f90"
echo "Source dir: $SRC"

cd "$SRC"
source ../../setup-env-channel.sh --cluster
makeneko src/turb_channel_two_phase.f90

mv neko neko_two_phase
echo "Done: $SRC/neko_two_phase"
