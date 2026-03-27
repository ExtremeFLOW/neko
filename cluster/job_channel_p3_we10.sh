#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 06:00:00
#SBATCH -N 2
#SBATCH -J channel_p3_we10
#SBATCH -p main
#SBATCH -o /cfs/klemming/projects/supr/kthmech/eriksie/logs/%J_channel_p3_we10.log

# Phase 3 two-phase channel: We=10 (σ=0.04), 144×18×48 mesh, t=20→25
#
# We=10: σ=0.04, Δt_cap≈0.079 (at Δx=0.0873), Δt/Δt_cap≈0.051
# Finest mesh in convergence study: 4.1 elements across interface.
#
# Prerequisites:
#   1. neko_two_phase built: sbatch build_neko_two_phase.sh
#   2. channel_p3_single_phase spin-up completed (fluid00004.chkp at t=20)
#   3. p3_sigma0 completed (confirm κ_rms level)
#
# Submit from Dardel login node:
#   sbatch $KTHMECH_PROJECT/scripts/job_channel_p3_we10.sh

set -e

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

SRC=$KTHMECH_PROJECT/src/neko-multiphase-channel/examples/two_phase_channel
RUN_NAME="channel_p3_we10"
RUN_DIR=$SCRATCH_DIR/$RUN_NAME

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

cp "$SRC/cases/144x24x48/turb_channel_two_phase_p3_we10.case" .
cp "$SRC/box_phys_144x24x48.nmsh" .
cp "$SRC/neko_two_phase" neko

# Copy checkpoint from spin-up if not already present
SPINUP_CHKP=$SCRATCH_DIR/channel_p3_single_phase/fluid00004.chkp
if [ ! -f fluid00004.chkp ]; then
    if [ -f "$SPINUP_CHKP" ]; then
        cp "$SPINUP_CHKP" .
        echo "Copied fluid00004.chkp from channel_p3_single_phase"
    else
        echo "ERROR: fluid00004.chkp not found. Run channel_p3_single_phase first."
        exit 1
    fi
fi

echo "Starting: $RUN_NAME"
echo "Run dir:  $RUN_DIR"
echo "Nodes: $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS"

srun -u -n 256 ./neko turb_channel_two_phase_p3_we10.case
