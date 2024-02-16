#!/usr/bin/bash

# ============================================================================ #
# Define the help function

function help() {
    echo -e "run.sh [Nx] [Ny] [Nz]"
    echo -e "  Generate a mesh and run the immersed_bunny case."
    echo -e "  The input arguments are the number of cells in the x, y, and z"
    echo -e "  directions, respectively."
    echo -e ""
    echo -e "  If no input arguments are provided, the default mesh size is"
    echo -e "  16x8x8."
    echo -e ""
    echo -e "  See Readme for additional details."
    exit 0
}

# Search for "-h" or "--help" in the input arguments
for arg in $@; do
    if [[ $arg == "-h" ]] || [[ $arg == "--help" ]]; then
        help
    fi
done

# ============================================================================ #
# Ensure Neko can be found and set default mesh size

if [ ! -z $NEKO_DIR ]; then
    export PATH=$NEKO_DIR/bin:$PATH
elif [[ -z $(which neko) ]]; then
    NEKO_DIR=$(realpath $0 | xargs dirname)
    NEKO_DIR=${NEKO_DIR%/example*}/external/neko
    export PATH=$NEKO_DIR/bin:$PATH
fi

if [ $# -lt 1 ]; then Nx=16; else Nx=$1; fi
if [ $# -lt 2 ]; then Ny=8; else Ny=$2; fi
if [ $# -lt 3 ]; then Nz=8; else Nz=$3; fi

# ============================================================================ #
# Generate mesh and run case

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox -100 400 -100 100 0 200 $Nx $Ny $Nz .false. .false. .false.
neko immersed_bunny.case

# End of file
# ============================================================================ #
