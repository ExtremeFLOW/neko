#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

"${PYTHON:-python3}" ../prepare_case.py --check
if [[ ! -x ./neko ]]; then
    printf '%s\n' 'Build first: makeneko sodermalm_wind.f90' >&2
    exit 1
fi
if [[ -d fields ]] && [[ -n $(find fields -type f -print -quit) ]]; then
    printf '%s\n' 'Refusing to overwrite existing fields; archive the previous run first.' >&2
    exit 1
fi

# Pass the site-specific launcher and its arguments, preserving GPU binding.
if [[ $# -eq 0 ]]; then
    printf '%s\n' 'Usage: bash run.sh LAUNCHER [ARGS...] (e.g. mpirun -np 6)' >&2
    exit 1
fi
exec "$@" ./neko sharp_mask.case
