#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -J build_neko_two_phase_p2
#SBATCH -p main
#SBATCH -o /cfs/klemming/projects/supr/kthmech/eriksie/logs/%J_build_neko_two_phase_p2.log

# Compile the _p2 two-phase user module on Dardel.
# This variant performs an extra GS pass on n̂ before div(n) to smooth
# the C0 kink at element faces (Lagrange endpoint derivative artifact).
#
# Run AFTER build_neko_channel.sh has completed (makeneko must be installed).
#
# Output: $SRC/neko_two_phase_p2  (used by job_channel_p3_sigma0_p2.sh etc.)
#
# Submit from the Dardel login node:
#   sbatch $KTHMECH_PROJECT/scripts/build_neko_two_phase_p2.sh

set -e

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

SRC=$KTHMECH_PROJECT/src/neko-multiphase-channel/examples/two_phase_channel

echo "Compiling turb_channel_two_phase_p2.f90"
echo "Source dir: $SRC"

cd "$SRC"
source ../../setup-env-channel.sh --cluster
makeneko src/turb_channel_two_phase_p2.f90

mv neko neko_two_phase_p2
echo "Done: $SRC/neko_two_phase_p2"
