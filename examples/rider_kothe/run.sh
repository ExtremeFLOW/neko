#!/bin/bash
# Usage: ./run.sh [case ...]
#
# Mesh, generated once -- not part of this script's automatic flow:
#   genmeshbox 0 1 0 1 -0.1 0 64 64 1 .true. .true. .true.
# 64x64 elements on a unit periodic square; with polynomial_order 5 that is
# H/N = 3.125e-3, so epsilon 3.125e-3 is xi = 1.0.
#
# dt is set by the CDI compression CFL gamma*u_max*dt/h_gll_min <= 0.05, not by
# the advective CFL Saini's own dt = 4e-4 assumes -- the compression limit is
# the tighter of the two here by more than an order of magnitude.
#
# The full run is t = 8, i.e. 100000 steps: launch detached.
#   setsid nohup ./run.sh > chain.log 2>&1 < /dev/null &
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -f box.nmsh ]; then
  echo "box.nmsh not found -- generate it first, see the comment above." >&2
  exit 1
fi

makeneko rider_kothe.f90

cases=("$@")
if [ ${#cases[@]} -eq 0 ]; then
  # the recommended setting first; the upstream baseline is kept so the one
  # prior run remains reproducible
  cases=(rider_kothe rider_kothe_xi10)
fi

if [ -n "${NEKO_PSI_TRANSPORT_CUDA_PREFIX:-}" ] &&
   [ "$(command -v neko)" = "$NEKO_PSI_TRANSPORT_CUDA_PREFIX/bin/neko" ]; then
  for c in "${cases[@]}"; do
    c="${c%.case}"
    echo "=== ${c} (gpu) ==="
    mpirun -np 1 ./neko "${c}.case" > "run_${c}.log" 2>&1
  done
else
  for c in "${cases[@]}"; do
    c="${c%.case}"
    echo "=== ${c} (cpu) ==="
    mpirun -np 4 --bind-to none ./neko "${c}.case" > "run_${c}.log" 2>&1
  done
fi
echo done
