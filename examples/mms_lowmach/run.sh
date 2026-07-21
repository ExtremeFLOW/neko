#!/usr/bin/env bash
# Build and run the low-Mach momentum-stress MMS verification.
#
# Verifies lowmach_vel_res_cpu_compute by feeding it analytic fields whose
# divergence of the deviatoric stress is known in closed form, then checking
# spectral convergence of the discrete result under p-refinement.
set -euo pipefail

NEKO=/home/hochi/neko-install
SRC=/home/hochi/neko-lowmach/src
JF=/home/hochi/json-fortran-install
export LD_LIBRARY_PATH=$JF/lib

cd "$(dirname "$0")"

# 1. Fully periodic box [0,2pi]^3 with 4^3 elements.
"$NEKO/bin/genmeshbox" \
    0 6.283185307179586 0 6.283185307179586 0 6.283185307179586 \
    4 4 4 .true. .true. .true.

# 2. Compile the standalone driver against the built libneko.
mpif90 -g -O2 -I"$SRC" mms_lowmach.f90 "$SRC/.libs/libneko.a" \
    -L"$JF/lib" -ljsonfortran -lopenblas -o mms_lowmach

# 3. Run (single rank).
./mms_lowmach
