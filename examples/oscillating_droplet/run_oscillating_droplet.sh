#!/bin/bash
# Deploy and run the oscillating droplet case on Dardel.
#
# Usage:
#   ./run_oscillating_droplet.sh              # local: run neko directly
#   ./run_oscillating_droplet.sh --cluster    # Dardel: submit via sbatch
#
# Prerequisites:
#   - box.nmsh in this directory (40x40x1 mesh, same as static droplet)
#   - oscillating_droplet.case and oscillating_droplet.f90 in this directory
#   - Cluster: dardel_job.sh template in this directory; sbatch available
#
# Output:
#   - Local:   ekin.csv + field files in this directory
#   - Cluster: /cfs/klemming/scratch/e/eriksie/oscillating_droplet/

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRATCH_DIR=/cfs/klemming/scratch/e/eriksie

# Parse arguments
CLUSTER_MODE=false
if [ "$1" = "--cluster" ]; then
    CLUSTER_MODE=true
elif [ -n "$1" ]; then
    NEKO="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
elif [ -f "$SCRIPT_DIR/neko" ]; then
    NEKO="$SCRIPT_DIR/neko"
else
    NEKO="$(command -v neko)"
fi

if $CLUSTER_MODE; then
    echo "Mode: cluster (sbatch)"
    if [ ! -f "$SCRIPT_DIR/dardel_job.sh" ]; then
        echo "ERROR: dardel_job.sh not found in $SCRIPT_DIR"
        exit 1
    fi
    RUN_DIR="$SCRATCH_DIR/oscillating_droplet"
    echo "Output: $RUN_DIR"

    mkdir -p "$RUN_DIR"

    # Copy required files to scratch
    cp -f "$SCRIPT_DIR/box.nmsh" "$RUN_DIR/"
    cp -f "$SCRIPT_DIR/oscillating_droplet.case" "$RUN_DIR/"
    cp -f "$SCRIPT_DIR/oscillating_droplet.f90" "$RUN_DIR/"
    cp -f "$SCRIPT_DIR/dardel_job.sh" "$RUN_DIR/job.sh"

    # Submit
    job_id=$(cd "$RUN_DIR" && sbatch --parsable job.sh)
    echo "Submitted job $job_id"
    echo "Monitor with: squeue -u \$USER"
    echo "Results in: $RUN_DIR"
else
    echo "Mode: local  (neko: $NEKO)"
    (cd "$SCRIPT_DIR" && "$NEKO" oscillating_droplet.case)
    echo "Done."
fi
