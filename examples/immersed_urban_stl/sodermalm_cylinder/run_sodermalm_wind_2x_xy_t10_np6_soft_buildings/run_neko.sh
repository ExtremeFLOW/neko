#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
export DYLD_LIBRARY_PATH="/usr/local/neko_install/bin:/usr/local/json-fortran_install/lib:${DYLD_LIBRARY_PATH:-}"

NP="${NP:-6}"
MPIEXEC="${MPIEXEC:-/opt/homebrew/bin/mpirun}"
NEKO_BIN="${NEKO_BIN:-/Users/tadao/KTH/neko-2/.codex_build/cpu_install/bin/neko}"
MAKENEKO="${MAKENEKO:-/Users/tadao/KTH/neko-2/.codex_build/cpu_install/bin/makeneko}"
PYTHON="${PYTHON:-/usr/local/bin/python3}"
CASE="sodermalm_wind_2x_xy_t10_np6_soft_buildings.case"
GENERATED="../generated_p7_xy35_inlet220_nz10"
MESH="$GENERATED/sodermalm_cylinder_p7_xy35_inlet220_nz10.nmsh"
USER_FILE="sodermalm_wind_profile_prefilled.f90"
CACHE_DIR="../cache_p3_xy35_inlet220_nz10_allbuildings_centered_sw50_h65"
CACHE_PREFIX="$CACHE_DIR/bldp3_ctr50"
CACHE_DATA="${CACHE_PREFIX}0.f00000"
CACHE_INDEX="${CACHE_PREFIX}0.nek5000"
TEMPLATE_FIELD="$CACHE_DIR/template.f00000"

export SODERMALM_INITIAL_WIND_SCALE="${SODERMALM_INITIAL_WIND_SCALE:-1.0}"
export SODERMALM_RAMP_TIME="${SODERMALM_RAMP_TIME:-0.0}"
export SODERMALM_INLET_ARC_WIDTH_DEG="${SODERMALM_INLET_ARC_WIDTH_DEG:-220.0}"
export SODERMALM_INLET_TAPER_DEG="${SODERMALM_INLET_TAPER_DEG:-0.0}"

if [[ ! -f "$MESH" ]]; then
  echo "Missing mesh: $MESH"
  echo "Generate the p7_xy35_inlet220_nz10 mesh first; see ../README.md."
  exit 1
fi

mkdir -p fields renders

if [[ ! -f "$CACHE_DATA" || ! -f "$CACHE_INDEX" ]]; then
  mkdir -p "$CACHE_DIR"
  "$PYTHON" ../make_template_field_from_nmsh.py \
    "$MESH" "$TEMPLATE_FIELD" --polynomial-order 3
  "$PYTHON" ../make_sodermalm_building_cache.py \
    --generated "$GENERATED" \
    --template-field "$TEMPLATE_FIELD" \
    --cache-prefix "$CACHE_PREFIX" \
    --smooth-width 50 \
    --transition-placement centered \
    --max-height 65
fi

cp "$CACHE_DATA" fields/
cp "$CACHE_INDEX" fields/

if [[ -x "$MAKENEKO" ]]; then
  "$MAKENEKO" "$USER_FILE"
  NEKO_RUN_BIN="./neko"
else
  echo "Warning: makeneko not found at $MAKENEKO; using NEKO_BIN=$NEKO_BIN"
  NEKO_RUN_BIN="$NEKO_BIN"
fi

{
  echo "Starting qualified Södermalm p3 diagnostic at $(date)"
  echo "Run directory: $(pwd)"
  echo "Case: $CASE"
  echo "Ranks: $NP"
  "$MPIEXEC" -np "$NP" "$NEKO_RUN_BIN" "$CASE"
  echo "Finished at $(date)"
} 2>&1 | tee run.log
