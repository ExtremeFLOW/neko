#!/bin/bash
#SBATCH -A naiss2025-3-39
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -J build_neko_channel
#SBATCH -p main
#SBATCH -o %J_build_neko_channel.log

# Build neko-multiphase-channel on Dardel.
# Submit from $KTHMECH_PROJECT/scripts/:
#   sbatch $KTHMECH_PROJECT/scripts/build_neko_channel.sh
#
# Source is expected at: $KTHMECH_PROJECT/src/neko-multiphase-channel/
# (clone with: cd $KTHMECH_PROJECT/src
#              git clone git@github.com:ExtremeFLOW/neko.git neko-multiphase-channel
#              cd neko-multiphase-channel && git checkout eriksie/multiphase/two-phase-channel)

set -e

SRC=$KTHMECH_PROJECT/src/neko-multiphase-channel

echo "Building from: $SRC"
echo "Installing to: $KTHMECH_PROJECT/builds/neko-channel"

cd "$SRC"
bash build-neko-channel.sh --cluster --clean
