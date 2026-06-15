#!/usr/bin/env bash
# set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
WORK_DIR="${1:-$SCRIPT_DIR/workdir}"
EXAMPLE_DIR="$REPO_ROOT/examples/tgv"

if ! command -v valgrind >/dev/null 2>&1; then
    echo "Skipping valgrind regression test: valgrind not found in PATH."
    exit 77
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "python3 is required for valgrind regression checks." >&2
    exit 2
fi

if ! command -v makeneko >/dev/null 2>&1; then
    echo "makeneko is required and must be available on PATH." >&2
    exit 2
fi

for file in \
    "$SCRIPT_DIR/valgrind_runtime.supp" \
    "$SCRIPT_DIR/check_valgrind.py" \
    "$EXAMPLE_DIR/tgv.f90" \
    "$EXAMPLE_DIR/tgv.case" \
    "$EXAMPLE_DIR/512.nmsh"; do
    if [[ ! -f "$file" ]]; then
        echo "Required file not found: $file" >&2
        exit 2
    fi
done

mkdir -p "$WORK_DIR"
cp -f "$SCRIPT_DIR/valgrind_runtime.supp" "$WORK_DIR/"
cp -f "$SCRIPT_DIR/check_valgrind.py" "$WORK_DIR/"
cp -f "$EXAMPLE_DIR/tgv.case" "$WORK_DIR/valgrind_tgv.case"
cp -f "$EXAMPLE_DIR/tgv.f90" "$WORK_DIR/valgrind_tgv.f90"
cp -f "$EXAMPLE_DIR/512.nmsh" "$WORK_DIR/"

pushd "$WORK_DIR" >/dev/null

rm -f neko user.mod user.smod valgrind_regression.log valgrind_stdout.log valgrind_stderr.log

# Change the end time to run 10 iterations
dt=$(grep -E '"timestep"\s*:' valgrind_tgv.case | sed -E 's/.*"timestep"\s*:\s*([0-9eE.+-]+).*/\1/')
if [[ -z "$dt" ]]; then
    echo "Failed to extract timestep from valgrind_tgv.case" >&2
    exit 2
fi
end_time=$(python3 -c "print(10 * $dt)")
sed -i.bak -E "s/(\"end_time\"\s*:\s*)[0-9eE.+-]+/\1$end_time/" valgrind_tgv.case

makeneko valgrind_tgv.f90 || {
    echo "makeneko failed to compile valgrind_tgv.f90; see above for errors." >&2
    exit 2
}

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

valgrind \
    --leak-check=full \
    --show-leak-kinds=all \
    --show-error-list=no \
    --undef-value-errors=no \
    --num-callers=20 \
    --suppressions=valgrind_runtime.supp \
    --log-file=valgrind_regression.log \
    ./neko valgrind_tgv.case > valgrind_stdout.log 2> valgrind_stderr.log
valgrind_exit=$?

if [[ $valgrind_exit -ne 0 ]]; then
    echo "Neko exited with status $valgrind_exit under Valgrind; check valgrind_stderr.log for details." >&2
    exit 2
fi

if [[ ! -s valgrind_regression.log ]]; then
    echo "Valgrind did not produce a log file; check valgrind_stderr.log for details." >&2
    exit 2
fi

backend_type="$(grep -E 'Bcknd type:' valgrind_stdout.log | tail -n 1 | awk -F ':' '{print toupper($2)}' | xargs)"

# Use backend-specific defaults inferred from current regression baselines.
# Users can always override these through NEKO_VALGRIND_MAX_* environment vars.
case "$backend_type" in
    "CPU")
        default_definite="${NEKO_VALGRIND_MAX_DEFINITE_CPU:-0}"
        default_indirect="${NEKO_VALGRIND_MAX_INDIRECT_CPU:-0}"
        default_possible="${NEKO_VALGRIND_MAX_POSSIBLE_CPU:-0}"
        default_reachable="${NEKO_VALGRIND_MAX_REACHABLE_CPU:-1200000}"
        ;;
    "CPU (LIBXSMM)")
        default_definite="${NEKO_VALGRIND_MAX_DEFINITE_LIBXSMM:-0}"
        default_indirect="${NEKO_VALGRIND_MAX_INDIRECT_LIBXSMM:-0}"
        default_possible="${NEKO_VALGRIND_MAX_POSSIBLE_LIBXSMM:-0}"
        default_reachable="${NEKO_VALGRIND_MAX_REACHABLE_LIBXSMM:-1200000}"
        ;;
    "SX-AURORA")
        default_definite="${NEKO_VALGRIND_MAX_DEFINITE_SX:-0}"
        default_indirect="${NEKO_VALGRIND_MAX_INDIRECT_SX:-0}"
        default_possible="${NEKO_VALGRIND_MAX_POSSIBLE_SX:-0}"
        default_reachable="${NEKO_VALGRIND_MAX_REACHABLE_SX:-1200000}"
        ;;
    "ACCELERATOR (CUDA)")
        default_definite="${NEKO_VALGRIND_MAX_DEFINITE_CUDA:-0}"
        default_indirect="${NEKO_VALGRIND_MAX_INDIRECT_CUDA:-0}"
        default_possible="${NEKO_VALGRIND_MAX_POSSIBLE_CUDA:-3000}"
        default_reachable="${NEKO_VALGRIND_MAX_REACHABLE_CUDA:-9000000}"
        ;;
    "ACCELERATOR (HIP)")
        default_definite="${NEKO_VALGRIND_MAX_DEFINITE_HIP:-0}"
        default_indirect="${NEKO_VALGRIND_MAX_INDIRECT_HIP:-0}"
        default_possible="${NEKO_VALGRIND_MAX_POSSIBLE_HIP:-3000}"
        default_reachable="${NEKO_VALGRIND_MAX_REACHABLE_HIP:-9000000}"
        ;;
    "ACCELERATOR (OPENCL)")
        default_definite="${NEKO_VALGRIND_MAX_DEFINITE_OPENCL:-0}"
        default_indirect="${NEKO_VALGRIND_MAX_INDIRECT_OPENCL:-0}"
        default_possible="${NEKO_VALGRIND_MAX_POSSIBLE_OPENCL:-3000}"
        default_reachable="${NEKO_VALGRIND_MAX_REACHABLE_OPENCL:-9000000}"
        ;;
    *)
        default_definite="${NEKO_VALGRIND_MAX_DEFINITE_UNKNOWN:-0}"
        default_indirect="${NEKO_VALGRIND_MAX_INDIRECT_UNKNOWN:-0}"
        default_possible="${NEKO_VALGRIND_MAX_POSSIBLE_UNKNOWN:-0}"
        default_reachable="${NEKO_VALGRIND_MAX_REACHABLE_UNKNOWN:-1200000}"
        ;;
esac

if [[ -z "$backend_type" ]]; then
    backend_type="UNKNOWN"
fi

max_definite="${NEKO_VALGRIND_MAX_DEFINITE:-$default_definite}"
max_indirect="${NEKO_VALGRIND_MAX_INDIRECT:-$default_indirect}"
max_possible="${NEKO_VALGRIND_MAX_POSSIBLE:-$default_possible}"
max_reachable="${NEKO_VALGRIND_MAX_REACHABLE:-$default_reachable}"

echo "Valgrind limits for backend '$backend_type'"
echo "  definitely lost: $max_definite bytes"
echo "  indirectly lost: $max_indirect bytes"
echo "  possibly lost:   $max_possible bytes"
echo "  still reachable: $max_reachable bytes"

python3 ./check_valgrind.py \
    --log valgrind_regression.log \
    --max-definite-bytes "$max_definite" \
    --max-indirect-bytes "$max_indirect" \
    --max-possible-bytes "$max_possible" \
    --max-reachable-bytes "$max_reachable"

if [[ $? -ne 0 ]]; then
    echo "Valgrind regression test failed; see above for details." >&2
    exit 2
fi

echo "Valgrind regression test passed."
popd >/dev/null
