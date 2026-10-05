#!/bin/bash
# Usage: ./run.sh [case ...]   (no argument runs all four)
#
# Mesh, generated once -- not part of this script's automatic flow:
#   genmeshbox 0 1 0 1 -0.1 0 50 50 1 .true. .true. .true.
# 50x50 elements on a unit periodic square; with polynomial_order 5 that is
# H/N = 0.004, so epsilon 0.0112 is xi = 2.8.
#
# The rotation is u = pi(0.5-y), v = pi(x-0.5), i.e. one turn per t = 2, so
# end_time 20 is ten full rotations.
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

# Whichever setup-env*.sh you sourced decides the backend; makeneko links
# straight to ./neko and fails with "Text file busy" if a previous run still
# holds it, so rename it first rather than waiting.
makeneko zalesak_disk.f90

cases=("$@")
if [ ${#cases[@]} -eq 0 ]; then
  # the settled set. The resolution study (zalesak_p7, zalesak_p7_phi,
  # zalesak_h100, zalesak_h150 -- one rotation each at fixed xi=2.8) is run
  # by name, not by default: together they are ~5 h on a GPU.
  cases=(zalesak zalesak_phi_normal zalesak_redistance_phi zalesak_redistance_psi)
fi

# One GPU on this box, so a CUDA build takes the cases one at a time (~1.2 h
# each). The CPU build is ~10.6 h per case, so there it is worth giving each
# four ranks and running two at a time on sixteen cores.
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
