#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 08:00:00
#SBATCH -N 1
#SBATCH -J channel_p3_sp
#SBATCH -p main
#SBATCH -o /cfs/klemming/projects/supr/kthmech/eriksie/logs/%J_channel_p3_sp.log

# Phase 3 single-phase spin-up: 144x18x48 mesh, t=0->25, 128 ranks
#
# Mesh: 144x18x48 = 124416 elements (1.78x Phase 2), Delta_xz=0.0873
# 1 node / 128 ranks: 972 elem/rank
# Memory: ~10 GB (fits easily in 256 GB/node)
# Runtime: ~8.5h estimated (2.37x Phase 2's 3.6h)
#
# Purpose: establish turbulence on the Phase 3 mesh and produce fluid00004.chkp
# at t=20 for use as IC in all Phase 3 two-phase restart cases.
#
# Prerequisite: box_phys_144x18x48.nmsh must exist in SRC.
# Generate it on Dardel login node (not compute-intensive):
#   cd $SRC && genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 144 18 48 .true. .false. .true.
#   mv box.nmsh box_phys_144x18x48.nmsh
#
# Expected outputs in $SCRATCH_DIR/channel_p3_single_phase/:
#   ekin.csv            -- Ekin and u_max every 50 steps
#   fluid0000{0-5}.chkp -- checkpoints at t=0,5,10,15,20,25
#   field0.f0000{0-5}   -- velocity snapshots at t=5,10,15,20,25 + t=0
#   fluid00004.chkp     -- THE checkpoint for Phase 3 two-phase restarts (t=20)

set -e

module load PrgEnv-cray
module use /cfs/klemming/pdc/projects/hpcrd/modules
module load hpcrd
module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

SRC=/cfs/klemming/projects/supr/kthmech/eriksie/src/neko-multiphase-channel/examples/two_phase_channel
RUN=/cfs/klemming/scratch/e/eriksie/channel_p3_single_phase

mkdir -p $RUN
cd $RUN

cp $SRC/turb_channel_single_phase_p3.case .
cp $SRC/box_phys_144x18x48.nmsh .
cp $SRC/neko_single_phase neko

echo "Starting Phase 3 single-phase spin-up"
echo "Run dir: $RUN"
echo "Nodes: $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS"
echo "Case: turb_channel_single_phase_p3.case"
echo "Mesh: box_phys_144x18x48.nmsh (144x18x48 = 124416 elements)"

srun -u -n 128 ./neko turb_channel_single_phase_p3.case
