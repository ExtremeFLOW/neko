#!/bin/bash
# Sync two-phase channel source files to Dardel.
#
# Usage:
#   ./cluster/sync_to_dardel.sh
#
# Transfers:
#   examples/two_phase_channel/    → Dardel source dir
#   src/ (neko library sources)    → Dardel source dir
#
# Prerequisites:
#   - SSH alias 'dardel-ftn' configured in ~/.ssh/config for the file transfer node
#   - Repo synced to $KTHMECH_PROJECT/src/neko-multiphase-channel/ on Dardel

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DARDEL_HOST="dardel-ftn"
DARDEL_SRC="/cfs/klemming/projects/supr/kthmech/eriksie/src/neko-multiphase-channel"

echo "=== Syncing to Dardel ==="
echo "Source:      $REPO_ROOT"
echo "Destination: ${DARDEL_HOST}:${DARDEL_SRC}"
echo ""

# Sync case files and user module
echo "--- examples/two_phase_channel/ ---"
rsync -avz --progress \
    --include="*/" \
    --include="*.f90" \
    --include="*.case" \
    --include="*.py" \
    --include="*.sh" \
    --exclude="__pycache__/" \
    --exclude="*.pyc" \
    --exclude="*.nmsh" \
    --exclude="*.chkp" \
    --exclude="field0.*" \
    --exclude="neko" \
    --exclude="ekin.csv" \
    --exclude="*.gif" \
    --exclude="*.png" \
    --exclude="*" \
    "$REPO_ROOT/examples/two_phase_channel/" \
    "${DARDEL_HOST}:${DARDEL_SRC}/examples/two_phase_channel/"

# Sync modified library sources (e.g. chkp_file.f90, case.f90)
echo "--- src/ (library patches) ---"
rsync -avz --progress \
    --include="*.f90" \
    --include="*.h" \
    --include="*/" \
    --exclude="*" \
    "$REPO_ROOT/src/" \
    "${DARDEL_HOST}:${DARDEL_SRC}/src/"

# Sync cluster scripts
echo "--- cluster/ ---"
rsync -avz --progress \
    "$REPO_ROOT/cluster/" \
    "${DARDEL_HOST}:${DARDEL_SRC}/cluster/"

echo ""
echo "Sync complete."
echo "To update \$KTHMECH_PROJECT/scripts/ on Dardel, run:"
echo "  ssh dardel 'cp \$KTHMECH_PROJECT/src/neko-multiphase-channel/cluster/job_*.sh \$KTHMECH_PROJECT/scripts/'"
