#!/usr/bin/env bash
# Build json-fortran + the low-Mach Neko fork on an HPC login/build node (CPU).
# Adjust NEKO_DIR / module loads for your cluster, then:  bash build_neko_hpc.sh
set -euo pipefail

NEKO_DIR="${NEKO_DIR:-$HOME/neko-lowmach}"
PREFIX="$NEKO_DIR/install"
JSONFORTRAN="${JSONFORTRAN:-$HOME/json-fortran-install}"

# --- load your toolchain (EDIT for your site) ---------------------------------
# module load gcc openmpi openblas cmake        # example
: "${FC:=mpif90}"
: "${CC:=mpicc}"
echo "Using FC=$FC  CC=$CC"

# --- 1. json-fortran (small CMake project) ------------------------------------
if [ ! -f "$JSONFORTRAN/lib/libjsonfortran.so" ] && [ ! -f "$JSONFORTRAN/lib/libjsonfortran.a" ]; then
  echo "=== building json-fortran -> $JSONFORTRAN ==="
  tmp=$(mktemp -d)
  git clone --depth 1 https://github.com/jacobwilliams/json-fortran "$tmp/jf"
  cmake -S "$tmp/jf" -B "$tmp/jf/build" -DCMAKE_INSTALL_PREFIX="$JSONFORTRAN" \
        -DCMAKE_Fortran_COMPILER=gfortran -DUSE_GNU_INSTALL_CONVENTION=ON \
        -DJSON_FORTRAN_USE_THREADS=OFF -DENABLE_TESTS=OFF
  cmake --build "$tmp/jf/build" --target install -j
  rm -rf "$tmp"
fi
export PKG_CONFIG_PATH="$JSONFORTRAN/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="$JSONFORTRAN/lib:${LD_LIBRARY_PATH:-}"

# --- 2. Neko (autotools; CPU only -- do NOT enable cuda/hip) -------------------
echo "=== configuring + building Neko -> $PREFIX ==="
cd "$NEKO_DIR"
[ -x ./configure ] || ./regen.sh
./configure --prefix="$PREFIX" FC="$FC" CC="$CC" \
    PKG_CONFIG_PATH="$JSONFORTRAN/lib/pkgconfig"
make -j"$(nproc)"
make install

echo
echo "DONE. Build a case driver with:"
echo "  export LD_LIBRARY_PATH=$JSONFORTRAN/lib:\$LD_LIBRARY_PATH"
echo "  $PREFIX/bin/makeneko heated_hump.f90"
