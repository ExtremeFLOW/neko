#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
export DYLD_LIBRARY_PATH="/usr/local/neko_install/bin:/usr/local/json-fortran_install/lib:${DYLD_LIBRARY_PATH:-}"

NP="${NP:-6}"
MPIEXEC="${MPIEXEC:-/opt/homebrew/bin/mpirun}"
NEKO_BIN="${NEKO_BIN:-/Users/tadao/KTH/neko-2/.codex_build/cpu_install/bin/neko}"
CASE="sodermalm_wind_2x_xy_t10_np6_soft_buildings.case"
MESH="../generated_2x_xy/sodermalm_cylinder_2x_xy.nmsh"

if [[ ! -f "$MESH" ]]; then
  echo "Missing mesh: $MESH"
  echo "Run ../build_sodermalm_cylinder_2x_xy.py first."
  exit 1
fi

mkdir -p fields renders

{
  echo "Starting 2x-xy soft-buildings Södermalm wind run at $(date)"
  echo "Run directory: $(pwd)"
  echo "Case: $CASE"
  echo "Ranks: $NP"
  "$MPIEXEC" -np "$NP" "$NEKO_BIN" "$CASE"
  echo "Finished at $(date)"
} 2>&1 | tee run.log
