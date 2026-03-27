#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -J channel_p2_sigma0
#SBATCH -p main
#SBATCH -o /cfs/klemming/projects/supr/kthmech/eriksie/logs/%J_channel_p2_sigma0.log

# Phase 2 two-phase channel: σ=0 CDI diagnostic, 108×18×36 mesh, t=20→21
#
# Purpose: verify κ_rms on the finer mesh. Compare with Phase 1 baseline (κ_rms≈64).
# If coverage helps, expect lower κ_rms. Run this FIRST before we10/we1.
#
# Prerequisites:
#   1. neko_two_phase built: sbatch build_neko_two_phase.sh
#   2. channel_p2_single_phase spin-up completed (fluid00004.chkp at t=20)
#
# Submit from Dardel login node:
#   sbatch $KTHMECH_PROJECT/scripts/job_channel_p2_sigma0.sh

set -e

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

SRC=$KTHMECH_PROJECT/src/neko-multiphase-channel/examples/two_phase_channel
RUN_NAME="channel_p2_sigma0"
RUN_DIR=$SCRATCH_DIR/$RUN_NAME

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

cp "$SRC/cases/108x18x36/turb_channel_two_phase_p2_sigma0.case" .
cp "$SRC/box_phys_108x18x36.nmsh" .
cp "$SRC/neko_two_phase" neko

# Copy checkpoint from spin-up if not already present
SPINUP_CHKP=$SCRATCH_DIR/channel_p2_single_phase/fluid00004.chkp
if [ ! -f fluid00004.chkp ]; then
    if [ -f "$SPINUP_CHKP" ]; then
        cp "$SPINUP_CHKP" .
        echo "Copied fluid00004.chkp from channel_p2_single_phase"
    else
        echo "ERROR: fluid00004.chkp not found. Run channel_p2_single_phase first."
        exit 1
    fi
fi

echo "Starting: $RUN_NAME"
echo "Run dir:  $RUN_DIR"
echo "Nodes: $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS"

srun -u -n 128 ./neko turb_channel_two_phase_p2_sigma0.case
