#!/bin/bash
# Environment for building neko-multiphase-psi-transport on this workstation (egidius).
#
# Usage: source setup-env.sh
#
# Reuses the JSON-Fortran install already built for the other
# neko-multiphase-* sandboxes instead of rebuilding it.

export JSON_INSTALL="/lscratch/sieburgh/local/jsonfortran-gnu-9.2.1"
export NEKO_PSI_TRANSPORT_PREFIX="/lscratch/sieburgh/local/neko-multiphase-psi-transport"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}${JSON_INSTALL}/lib/"
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH:+$PKG_CONFIG_PATH:}${JSON_INSTALL}/lib/pkgconfig"

export FC="gfortran"
export CC="gcc"
export CXX="g++"

export PATH="$NEKO_PSI_TRANSPORT_PREFIX/bin:$PATH"

echo "neko-multiphase-psi-transport environment loaded (prefix: $NEKO_PSI_TRANSPORT_PREFIX)"
