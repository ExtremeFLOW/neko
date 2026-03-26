#!/bin/bash
# Sync simulation results from Dardel scratch to egidius.
#
# Usage:
#   ./cluster/sync_from_dardel.sh <run_name> [--all]
#
#   <run_name>  : directory under $SCRATCH_DIR on Dardel (e.g. channel_p2_sigma0)
#   --all       : sync ALL field files (default: only ekin.csv + field0.f00000 + sparse snapshots)
#
# Results land in: /lscratch/sieburgh/simulations/<run_name>/

set -e

if [ -z "$1" ]; then
    echo "Usage: $0 <run_name> [--all]"
    exit 1
fi

RUN_NAME="$1"
SYNC_ALL=false
if [ "$2" = "--all" ]; then
    SYNC_ALL=true
fi

DARDEL_HOST="dardel-ftn"
DARDEL_SCRATCH="/cfs/klemming/scratch/e/eriksie"
LOCAL_SIMS="/lscratch/sieburgh/simulations"

SRC="${DARDEL_HOST}:${DARDEL_SCRATCH}/${RUN_NAME}/"
DST="${LOCAL_SIMS}/${RUN_NAME}/"

mkdir -p "$DST"

echo "=== Syncing from Dardel ==="
echo "Source:      $SRC"
echo "Destination: $DST"
echo ""

if $SYNC_ALL; then
    echo "Mode: ALL files"
    rsync -avz --progress \
        --exclude="*.chkp" \
        --exclude="neko" \
        "$SRC" "$DST"
else
    echo "Mode: selective (ekin.csv + mesh + sparse snapshots)"

    # Always get the diagnostic CSV
    rsync -avz --progress \
        --include="ekin.csv" \
        --include="*.log" \
        --include="*.case" \
        --exclude="*" \
        "$SRC" "$DST"

    # Initial snapshot (mesh coordinates)
    rsync -avz --progress \
        --include="field0.f00000" \
        --exclude="field0.f*" \
        --exclude="*" \
        "$SRC" "$DST"

    # One snapshot per 0.5 TU (files ending in 0 or 5 in the last two digits)
    rsync -avz --progress \
        --include="field0.f*[05]" \
        --exclude="field0.f*" \
        --exclude="*" \
        "$SRC" "$DST"
fi

echo ""
echo "Done. Files in: $DST"
ls -lh "$DST"
