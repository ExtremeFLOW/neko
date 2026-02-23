#!/bin/bash
# Set up a local Python venv for spurious currents postprocessing.
# Uses python3.12 (required: pysemtools needs Python <3.13).
#
# Usage:
#   ./setup-python-env.sh          # create and populate venv (first time)
#   source .venv/bin/activate      # activate in future sessions

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$SCRIPT_DIR/.venv"
PYTHON=/usr/local/bin/python3.12
PYSEMTOOLS=/Volumes/sourcecode/Python/kth/pySEMTools

if [ ! -f "$PYTHON" ]; then
    echo "ERROR: python3.12 not found at $PYTHON"
    exit 1
fi

if [ ! -d "$PYSEMTOOLS" ]; then
    echo "ERROR: pysemtools source not found at $PYSEMTOOLS"
    exit 1
fi

echo "Creating venv with $($PYTHON --version) ..."
$PYTHON -m venv "$VENV"

echo "Installing pysemtools and dependencies ..."
"$VENV/bin/pip" install --quiet -e "$PYSEMTOOLS"

echo "Installing notebook packages ..."
"$VENV/bin/pip" install --quiet matplotlib jupyter ipykernel

echo ""
echo "Done. To activate:"
echo "  source .venv/bin/activate"
