#!/bin/bash
# Run spurious currents benchmark for multiple sigma values.
# Each run goes to t*=250, i.e., end_time = 10.0/sigma.
# Creates one subdirectory per sigma: sigma_10.0/, sigma_1.0/, etc.
#
# Usage:
#   ./run_sigma_sweep.sh              # local: run neko directly
#   ./run_sigma_sweep.sh --cluster    # Dardel: submit via sbatch
#   ./run_sigma_sweep.sh /path/neko   # local: explicit neko binary
#
# Prerequisites:
#   - box.nmsh in this directory (40x40x1 mesh)
#   - spurious_currents.case (template) in this directory
#   - Local: neko binary (makeneko-compiled) in this directory or on PATH
#   - Cluster: dardel_job.sh template in this directory; sbatch available
#
# Output:
#   - Local:   sigma_*/  subdirectories in this directory
#   - Cluster: /cfs/klemming/scratch/e/eriksie/spurious_currents_sigma_sweep/sigma_*/

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRATCH_DIR=/cfs/klemming/scratch/e/eriksie

# Parse arguments
CLUSTER_MODE=false
if [ "$1" = "--cluster" ]; then
    CLUSTER_MODE=true
elif [ -n "$1" ]; then
    # Explicit neko path — resolve to absolute so it stays valid after cd
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
    SWEEP_DIR="$SCRATCH_DIR/spurious_currents_sigma_sweep"
    echo "Each job will compile spurious_currents.f90 via makeneko"
    echo "Output: $SWEEP_DIR"
else
    SWEEP_DIR="$SCRIPT_DIR"
    echo "Mode: local  (neko: $NEKO)"
fi

SIGMA_VALUES=(10.0 1.0 0.5 0.1 0.05 0.01)

for sigma in "${SIGMA_VALUES[@]}"; do
    # Physical end time for t*=250: t = mu*D*250/sigma = 0.1*0.4*250/sigma = 10.0/sigma
    end_time=$(python3 -c "print(10.0 / $sigma)")

    dir="$SWEEP_DIR/sigma_${sigma}"
    echo "=== sigma=$sigma  end_time=$end_time  dir=$dir ==="

    mkdir -p "$dir"

    # Copy the mesh (symlinks may not work across filesystems)
    cp -f "$SCRIPT_DIR/box.nmsh" "$dir/box.nmsh"

    # Patch the case file: update sigma and end_time
    sed \
        -e "s/\"sigma\": [0-9.e+-]*/\"sigma\": $sigma/" \
        -e "s/\"end_time\": [0-9.e+-]*/\"end_time\": $end_time/" \
        "$SCRIPT_DIR/spurious_currents.case" > "$dir/spurious_currents.case"

    if $CLUSTER_MODE; then
        # Copy user file so each job compiles its own binary
        cp "$SCRIPT_DIR/spurious_currents.f90" "$dir/"

        # Instantiate the Slurm template and submit
        sed \
            -e "s/SIGMA/$sigma/g" \
            -e "s/END_TIME/$end_time/g" \
            "$SCRIPT_DIR/dardel_job.sh" > "$dir/job.sh"
        job_id=$(cd "$dir" && sbatch --parsable job.sh)
        echo "  Submitted job $job_id"
    else
        # Symlink the mesh locally (same filesystem)
        ln -sf "$SCRIPT_DIR/box.nmsh" "$dir/box.nmsh"
        # Run neko directly (local)
        (cd "$dir" && "$NEKO" spurious_currents.case)
        echo "=== Done sigma=$sigma ==="
    fi
done

if $CLUSTER_MODE; then
    echo ""
    echo "All jobs submitted. Monitor with: squeue -u \$USER"
    echo "Results in: $SWEEP_DIR"
else
    echo "All runs complete."
fi
