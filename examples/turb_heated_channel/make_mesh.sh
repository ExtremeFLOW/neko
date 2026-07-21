#!/usr/bin/env bash
# Generate the wall-clustered turbulent-heated-channel mesh -> box.nmsh
#
#   Domain : x in [0, 2*pi] (streamwise, periodic)
#            y in [0, 2]     (wall-normal, no-slip walls -> zones 3=y0, 4=y2)
#            z in [0, pi]    (spanwise,   periodic)
#   Elements: 16 (x) x 16 (y) x 12 (z), polynomial order 7 (set in the .case)
#   Wall-normal y is CLUSTERED toward both walls via ydist.csv (tanh, gamma=2);
#   x and z stay uniform. genmeshbox itself makes only uniform boxes, but it
#   accepts a per-direction vertex-distribution file (args 13-15).
#
#   genmeshbox arg order:
#     x0 x1 y0 y1 z0 z1  nelx nely nelz  perx pery perz  distx disty distz
set -euo pipefail
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=/home/hochi/json-fortran-install/lib:${LD_LIBRARY_PATH:-}
GENMESHBOX=/home/hochi/neko-lowmach/install/bin/genmeshbox

rm -f genmeshbox*.log box.nmsh
"$GENMESHBOX" 0.0 6.283185307179586 0.0 2.0 0.0 3.141592653589793 \
    16 16 12 \
    .true. .false. .true. \
    uniform ydist.csv uniform

echo "Wrote $(pwd)/box.nmsh"
