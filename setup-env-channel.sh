#!/bin/bash
# Environment setup for neko-multiphase-channel
#
# Usage:
#   source setup-env-channel.sh --local     # macOS development
#   source setup-env-channel.sh --cluster   # HPC cluster (Dardel/Cray + ROCm)

if [ "$1" = "--local" ]; then
    # ── macOS local development ──────────────────────────────────────────
    export JSON_INSTALL="$HOME/local"
    export NEKO_CHANNEL_PREFIX="$HOME/local/neko-channel"

    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}${JSON_INSTALL}/lib/"
    export PKG_CONFIG_PATH="${PKG_CONFIG_PATH:+$PKG_CONFIG_PATH:}${JSON_INSTALL}/lib/pkgconfig"

    export FC="ccache gfortran"
    export CC="ccache clang"
    export CXX="ccache clang++"

    alias make='make -j$(sysctl -n hw.ncpu)'

    export PATH="$NEKO_CHANNEL_PREFIX/bin:$PATH"

    echo "neko-channel environment loaded (local/macOS)"

elif [ "$1" = "--cluster" ]; then
    # ── Dardel (Cray + ROCm/HIP) ────────────────────────────────────────
    export KTHMECH_PROJECT=/cfs/klemming/projects/supr/kthmech/eriksie
    export SCRATCH_DIR=/cfs/klemming/scratch/e/eriksie
    export NEKO_CHANNEL_PREFIX="$KTHMECH_PROJECT/builds/neko-channel"

    module load PrgEnv-cray
    module use /cfs/klemming/pdc/projects/hpcrd/modules
    module load hpcrd
    module load json-fortran/8.3.0-cce-18.0.1-bjoug3p

    export PATH="$NEKO_CHANNEL_PREFIX/bin:$PATH"

    echo "neko-channel environment loaded (cluster/Dardel)"

else
    echo "Usage: source setup-env-channel.sh --local|--cluster"
    return 1 2>/dev/null || exit 1
fi
