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

makeneko valgrind_tgv.f90

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

valgrind -s \
    --leak-check=full \
    --show-leak-kinds=all \
    --undef-value-errors=no \
    --num-callers=20 \
    --suppressions=valgrind_runtime.supp \
    --log-file=valgrind_regression.log \
    ./neko valgrind_tgv.case > valgrind_stdout.log 2> valgrind_stderr.log

python3 ./check_valgrind.py \
    --log valgrind_regression.log \
    --max-definite-bytes "${NEKO_VALGRIND_MAX_DEFINITE:-0}" \
    --max-indirect-bytes "${NEKO_VALGRIND_MAX_INDIRECT:-0}" \
    --max-possible-bytes "${NEKO_VALGRIND_MAX_POSSIBLE:-0}" \
    --max-reachable-bytes "${NEKO_VALGRIND_MAX_REACHABLE:-1300000}"

echo "Valgrind regression test passed."
popd >/dev/null
