#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -J channel_p2_CASE
#SBATCH -p main
#SBATCH -o %J_channel_p2_CASE.log

# Dardel job script for Phase 2 two-phase channel runs.
#
# Before submitting:
#   1. Replace CASE with the actual case name (e.g. sigma0, we10, we1)
#   2. Adjust -t (wall time) and -N (nodes) as needed
#   3. Ensure fluid00004.chkp exists in the run directory (from new-mesh spin-up)
#
# Submit from egidius:
#   ssh dardel "sbatch \$KTHMECH_PROJECT/scripts/job_channel_p2_CASE.sh"
#
# Single-phase spin-up (new mesh, run this first):
#   Use N=2, -t 08:00:00, srun -n 256, end_time=25.0
#   fluid00004.chkp is written at t=20 (checkpoint_value=5.0, first checkpoint)

set -e

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

SRC=$KTHMECH_PROJECT/src/neko-multiphase-channel/examples/two_phase_channel
RUN_NAME="channel_p2_CASE"
RUN_DIR=$SCRATCH_DIR/$RUN_NAME

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

# Copy inputs (makeneko must already have been run in SRC)
cp "$SRC/cases/108x18x36/turb_channel_two_phase_p2_CASE.case" .  # adjust mesh dir for p3
cp "$SRC/box_phys_108x18x36.nmsh" .
cp "$SRC/neko" .

# fluid00004.chkp must already be present (copied from spin-up run)
if [ ! -f fluid00004.chkp ]; then
    echo "ERROR: fluid00004.chkp not found in $RUN_DIR"
    echo "Run the new-mesh single-phase spin-up first and copy the checkpoint here."
    exit 1
fi

echo "Starting: $RUN_NAME"
echo "Case file: turb_channel_two_phase_p2_CASE.case"
echo "Nodes: $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS"

srun -u -n 128 ./neko turb_channel_two_phase_p2_CASE.case
