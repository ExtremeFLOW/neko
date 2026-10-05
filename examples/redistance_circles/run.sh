#!/bin/bash
# Usage: ./run.sh [case ...]   (no argument runs the four committed cells)
#
# Saini et al. (2026) section 4.4 -- standalone re-distancing around two
# intersecting circles. Each cell is one Eq. (44) relaxation to tau = 6; there
# is no flow and no timestepping, so the whole result is produced in the user
# initialize hook and printed to the log. The single output frame holds the
# relaxed field.
#
# Meshes, generated once -- not part of this script's automatic flow:
#   genmeshbox -2 2 -2 2 -0.1 0 20 20 1 .false. .false. .true.  && mv box.nmsh box20.nmsh
#   genmeshbox -2 2 -2 2 -0.1 0 40 40 1 .false. .false. .true.  && mv box.nmsh box40.nmsh
#   genmeshbox -2 2 -2 2 -0.1 0 80 80 1 .false. .false. .true.  && mv box.nmsh box80.nmsh
# Non-periodic in x and y is load-bearing: the initial condition is not
# periodic (phi_0 is +3.03 at x=-2 against +0.63 at x=+2 on y=0), and Eq. (44) would
# carry a seam jump inward at unit speed across the whole domain by tau = 6.
#
# The Fig. 12 grid beyond the committed cells (3 meshes x N = 3..8) is generated
# from the same six-field template rather than committed; the committed set is
# Table 2's four cells.
#
# For long chains, launch detached:
#   setsid nohup ./run.sh > chain.log 2>&1 < /dev/null &
set -uo pipefail   # NOT -e: a diverging cell is a result, not a script failure
cd "$(dirname "$0")"

[ -z "${NEKO_PSI_TRANSPORT_CUDA_PREFIX:-}" ] && . ../../setup-env-cuda.sh

for m in box20.nmsh box40.nmsh box80.nmsh; do
  [ -f "$m" ] || { echo "$m not found -- generate it first, see above." >&2; exit 1; }
done

# makeneko links straight to ./neko and fails with "Text file busy" if a
# previous run still holds it, so rename rather than wait. Moving the file is
# safe for a running process, which keeps its own inode.
[ -f neko ] && mv -f neko neko.prev
makeneko redistance_circles.f90 || exit 1

# GPU occupancy ONLY -- counting `neko` processes would also match
# advecting_slab_1d, which is a CPU build (N = 10 is above the CUDA lx <= 10 cap).
gpu_busy() {
  local n
  n=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -c .)
  [ "${n:-0}" -gt 0 ]
}
wait_for_gpu() {
  local w=0
  while gpu_busy; do
    [ $w -eq 0 ] && echo "    GPU busy -- waiting (60 s poll, giving up after 8 h)"
    sleep 60; w=$((w+60))
    if [ $w -ge 28800 ]; then echo "    !! GPU still busy after 8 h, skipping $1"; return 1; fi
  done
  [ $w -gt 0 ] && echo "    GPU free after ${w}s"
  return 0
}

cases=("$@")
if [ ${#cases[@]} -eq 0 ]; then
  # Table 2's four cells, cheapest first.
  cases=(circles_h10_n3 circles_h5_n7 circles_h20_n3 circles_h10_n7)
fi

for c in "${cases[@]}"; do
  c="${c%.case}"
  [ -f "$c.case" ] || { echo "!! no case file for $c" >&2; continue; }
  echo "=== $c  $(date +%H:%M:%S) ==="
  wait_for_gpu "$c" || continue
  t0=$SECONDS
  mpirun -np 1 ./neko "$c.case" > "run_$c.log" 2>&1
  rc=$?
  printf '    exit %d after %d s   %s\n' "$rc" "$((SECONDS-t0))" \
    "$(grep -a '/int psi_e:' "run_$c.log" | tail -1)"
done
echo done
