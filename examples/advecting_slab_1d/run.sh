#!/bin/bash
# Usage: ./run.sh [case ...]   (no argument runs the default seven below)
#
# Mesh, generated once -- not part of this script's automatic flow:
#   genmeshbox 0 1 0 0.1 0 0.1 10 1 1 .true. .true. .true.
# 10 elements across a unit-length periodic domain; with polynomial_order 10
# that is H/N = 0.01, so epsilon 0.01 is xi = 1.0 and 0.005 is xi = 0.5.
#
# For long runs, launch detached so closing the terminal does not kill them:
#   setsid nohup ./run.sh > chain.log 2>&1 < /dev/null &
# then poll run_<case>.log.
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -f box.nmsh ]; then
  echo "box.nmsh not found -- generate it first, see the comment above." >&2
  exit 1
fi

# makeneko links straight to ./neko and fails with "Text file busy" if a
# previous run still holds it. Rename it first; the running process keeps its
# own inode.
makeneko advecting_slab_1d.f90

cases=("$@")
if [ ${#cases[@]} -eq 0 ]; then
  # the four xi/normal runs, plus the gamma line at xi=1.0 that shows the
  # method beats being switched off (gamma=0 is the null: no compression,
  # no balancing diffusion, just advection)
  cases=(phi_xi10 psi_xi10 phi_xi05 psi_xi05 gamma0_xi10 gamma025_xi10 psi_xi28)
fi

# One GPU on this box, so a CUDA build runs the cases one at a time. The CPU
# build gives each case a rank and runs them together -- they are ~13 min each.
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
    mpirun -np 1 --bind-to none ./neko "${c}.case" > "run_${c}.log" 2>&1 &
  done
  wait
fi
echo done
