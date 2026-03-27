#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 16:00:00
#SBATCH -N 4
#SBATCH -J channel_l3_sp
#SBATCH -p main
#SBATCH -o /cfs/klemming/projects/supr/kthmech/eriksie/logs/%J_channel_l3_sp.log

# L3 single-phase spin-up: 192x32x64 mesh, t=0->25, 256 ranks
#
# Mesh: 192x32x64 = 393,216 elements, Δxz=0.0654
# 4 nodes / 512 ranks: 768 elem/rank (2 nodes/256 ranks OOM at 1536 elem/rank)
# Runtime: ~12h estimated
#
# Purpose: establish turbulence on the L3 mesh and produce fluid00004.chkp
# at t=20 for use as IC in all L3 two-phase restart cases.
#
# Prerequisite: box_phys_192x32x64.nmsh must exist in SRC (root of two_phase_channel/).
# Generate on egidius and transfer: rsync box_phys_192x32x64.nmsh dardel-ftn:$SRC/
#
# Expected outputs in $SCRATCH_DIR/channel_l3_single_phase/:
#   ekin.csv            -- Ekin and u_max every 50 steps
#   fluid0000{0-5}.chkp -- checkpoints at t=0,5,10,15,20,25
#   field0.f0000{0-5}   -- velocity snapshots at t=5,10,15,20,25 + t=0
#   fluid00004.chkp     -- THE checkpoint for L3 two-phase restarts (t=20)

set -e

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

SRC=/cfs/klemming/projects/supr/kthmech/eriksie/src/neko-multiphase-channel/examples/two_phase_channel
RUN=/cfs/klemming/scratch/e/eriksie/channel_l3_single_phase

mkdir -p $RUN
cd $RUN

cp $SRC/cases/192x32x64/turb_channel_single_phase_l3.case .
cp $SRC/box_phys_192x32x64.nmsh .
cp $SRC/neko_single_phase neko

echo "Starting L3 single-phase spin-up"
echo "Run dir: $RUN"
echo "Nodes: $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS"
echo "Case: turb_channel_single_phase_l3.case"
echo "Mesh: box_phys_192x32x64.nmsh (192x32x64 = 393216 elements)"

srun -u -n 512 ./neko turb_channel_single_phase_l3.case
