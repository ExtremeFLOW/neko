#!/bin/bash
# Environment for the CUDA (GPU) build of neko-multiphase-psi-transport on
# this workstation (egidius, NVIDIA RTX 3090, compute capability 8.6).
#
# Usage: source setup-env-cuda.sh
#
# Reuses the JSON-Fortran install already built for the other
# neko-multiphase-* sandboxes. CUDA itself is the system (Debian-packaged)
# toolkit under /usr, not a separate /usr/local/cuda tree.

export JSON_INSTALL="/lscratch/sieburgh/local/jsonfortran-gnu-9.2.1"
export NEKO_PSI_TRANSPORT_CUDA_PREFIX="/lscratch/sieburgh/local/neko-multiphase-psi-transport-cuda"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}${JSON_INSTALL}/lib/"
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH:+$PKG_CONFIG_PATH:}${JSON_INSTALL}/lib/pkgconfig"

export FC="gfortran"
export CC="gcc"
export CXX="g++"

export PATH="$NEKO_PSI_TRANSPORT_CUDA_PREFIX/bin:$PATH"

echo "neko-multiphase-psi-transport CUDA environment loaded (prefix: $NEKO_PSI_TRANSPORT_CUDA_PREFIX)"
