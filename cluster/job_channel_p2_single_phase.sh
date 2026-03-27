#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -J channel_p2_sp
#SBATCH -p main
#SBATCH -o /cfs/klemming/projects/supr/kthmech/eriksie/logs/%J_channel_p2_sp.log

# Phase 2 single-phase spin-up: 108x18x36 mesh, t=0->25, 128 ranks
#
# Purpose: establish turbulence on the new mesh and produce fluid00004.chkp
# at t=20 for use as IC in all Phase 2 two-phase restart cases.
#
# Expected outputs in $SCRATCH_DIR/channel_p2_single_phase/:
#   ekin.csv            -- Ekin and u_max every 50 steps (verify u_max ~ 1.35-1.45)
#   fluid0000{0-5}.chkp -- checkpoints at t=0,5,10,15,20,25
#   field0.f0000{0-5}   -- velocity snapshots at t=5,10,15,20,25 + t=0
#   fluid00004.chkp     -- THE checkpoint to use for two-phase restarts (t=20)

set -e

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

SRC=/cfs/klemming/projects/supr/kthmech/eriksie/src/neko-multiphase-channel/examples/two_phase_channel
RUN=/cfs/klemming/scratch/e/eriksie/channel_p2_single_phase

mkdir -p $RUN
cd $RUN

cp $SRC/turb_channel_single_phase_p2.case .
cp $SRC/box_phys_108x18x36.nmsh .
cp $SRC/neko_single_phase neko

echo "Starting single-phase spin-up"
echo "Run dir: $RUN"
echo "Nodes: $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS"
echo "Case: turb_channel_single_phase_p2.case"

srun -u -n 128 ./neko turb_channel_single_phase_p2.case
