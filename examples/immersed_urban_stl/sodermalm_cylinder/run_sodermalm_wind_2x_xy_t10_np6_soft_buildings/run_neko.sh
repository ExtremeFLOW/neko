#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
export DYLD_LIBRARY_PATH="/usr/local/neko_install/bin:/usr/local/json-fortran_install/lib:${DYLD_LIBRARY_PATH:-}"

NP="${NP:-6}"
MPIEXEC="${MPIEXEC:-/opt/homebrew/bin/mpirun}"
NEKO_BIN="${NEKO_BIN:-/Users/tadao/KTH/neko-2/.codex_build/cpu_install/bin/neko}"
MAKENEKO="${MAKENEKO:-/Users/tadao/KTH/neko-2/.codex_build/cpu_install/bin/makeneko}"
CASE="sodermalm_wind_2x_xy_t10_np6_soft_buildings.case"
MESH="../generated_2x_xy/sodermalm_cylinder_2x_xy.nmsh"
USER_FILE="sodermalm_wind_profile.f90"

if [[ ! -f "$MESH" ]]; then
  echo "Missing mesh: $MESH"
  echo "Run ../build_sodermalm_cylinder_2x_xy.py first."
  exit 1
fi

mkdir -p fields renders

if [[ -x "$MAKENEKO" ]]; then
  "$MAKENEKO" "$USER_FILE"
  NEKO_RUN_BIN="./neko"
else
  echo "Warning: makeneko not found at $MAKENEKO; using NEKO_BIN=$NEKO_BIN"
  NEKO_RUN_BIN="$NEKO_BIN"
fi

{
  echo "Starting 2x-xy soft-buildings Södermalm wind run at $(date)"
  echo "Run directory: $(pwd)"
  echo "Case: $CASE"
  echo "Ranks: $NP"
  "$MPIEXEC" -np "$NP" "$NEKO_RUN_BIN" "$CASE"
  echo "Finished at $(date)"
} 2>&1 | tee run.log
