#!/bin/bash
# Build script for neko-multiphase-channel
#
# Usage:
#   ./build-neko-channel.sh --local            # build for macOS
#   ./build-neko-channel.sh --local --clean     # full clean rebuild
#   ./build-neko-channel.sh --local --regen     # re-run regen.sh + reconfigure

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Determine environment ───────────────────────────────────────────────────
ENV_FLAG=""
BUILD_FLAG=""

for arg in "$@"; do
    case "$arg" in
        --local|--cluster) ENV_FLAG="$arg" ;;
        --clean|--regen)   BUILD_FLAG="$arg" ;;
    esac
done

if [ -z "$ENV_FLAG" ]; then
    echo "Usage: ./build-neko-channel.sh --local|--egidius|--cluster [--clean|--regen]"
    exit 1
fi

# Source the environment
source "$SCRIPT_DIR/setup-env-channel.sh" "$ENV_FLAG"

INSTALL_PREFIX="$NEKO_CHANNEL_PREFIX"

echo "=== Building neko-channel ==="
echo "Source:  $SCRIPT_DIR"
echo "Install: $INSTALL_PREFIX"
echo ""

cd "$SCRIPT_DIR"

# ── Generate build system ────────────────────────────────────────────────────
if [ ! -f configure ] || [ "$BUILD_FLAG" = "--regen" ] || [ "$BUILD_FLAG" = "--clean" ]; then
    echo "--- Running regen.sh ---"
    ./regen.sh
fi

# ── Clean if requested ───────────────────────────────────────────────────────
if [ "$BUILD_FLAG" = "--clean" ]; then
    echo "--- Cleaning previous build ---"
    make distclean 2>/dev/null || true
fi

# ── Configure ────────────────────────────────────────────────────────────────
if [ ! -f Makefile ] || [ "$BUILD_FLAG" = "--regen" ] || [ "$BUILD_FLAG" = "--clean" ]; then
    echo "--- Configuring ---"
    if [ "$ENV_FLAG" = "--local" ]; then
        ./configure --prefix="$INSTALL_PREFIX" \
            LDFLAGS="-L${JSON_INSTALL}/lib -Wl,-rpath,${JSON_INSTALL}/lib" \
            FCFLAGS="-O3 -march=native -g"
    elif [ "$ENV_FLAG" = "--egidius" ]; then
        ./configure --prefix="$INSTALL_PREFIX" \
            LDFLAGS="-L${JSON_INSTALL}/lib -Wl,-rpath,${JSON_INSTALL}/lib" \
            FCFLAGS="-O3 -march=native -g"
    elif [ "$ENV_FLAG" = "--cluster" ]; then
        ./configure --prefix="$INSTALL_PREFIX" \
            MPIFC=ftn MPICC=cc MPICXX=CC FC=ftn CC=cc \
            FCFLAGS="-m4 -O2"
    fi
fi

# ── Build and install ────────────────────────────────────────────────────────
echo "--- Building ---"
make -j 32 install

echo ""
echo "=== neko-channel installed to $INSTALL_PREFIX ==="
echo "Usage:"
echo "  $INSTALL_PREFIX/bin/makeneko your_case.f90"
echo "  mpirun -np 4 ./neko your_case.case"
