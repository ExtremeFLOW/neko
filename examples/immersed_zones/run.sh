#!/usr/bin/bash

# ============================================================================ #
# Define the help function

function help() {
    echo -e "run.sh case"
    echo -e "  Generate a mesh and run all the desired case."
    echo -e "  The input arguments are the number of cells in the x, y, and z"
    echo -e "  directions, respectively."
    echo -e ""
    echo -e "  If no input arguments are provided, the default mesh size is"
    echo -e "  16x8x8."
    echo -e ""
    echo -e "  Example usage:"
    echo -e "    run.sh 16 8 8"
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    echo -e "  -x          Number of cells in the x direction."
    echo -e "  -y          Number of cells in the y direction."
    echo -e "  -z          Number of cells in the z direction."
    echo -e ""
    echo -e "  See Readme for additional details."
    exit 0
}

Nx=16
Ny=8
Nz=8

# Search for "-h" or "--help" in the input arguments
for arg in $@; do
    case $arg in
    -h | --help) help ;;
    -x) Nx=$arg ;;
    -y) Ny=$arg ;;
    -z) Nz=$arg ;;
    *) ;;
    esac
done

# Extract the cases to run
cases=""
for arg in $@; do
    case $arg in
    -h | --help | -x | -y | -z) continue ;;
    *) cases+="$arg " ;;
    esac
done

# ============================================================================ #
# Ensure Neko can be found and set default mesh size

if [ "$NEKO_DIR" ]; then
    export PATH=$NEKO_DIR:$PATH
fi

if [[ -z $(which neko) ]]; then
    echo -e "Neko not found." >&2
    echo -e "Please ensure Neko is installed and in your PATH." >&2
    echo -e "Alternatively, set the NEKO_DIR environment variable." >&2
    exit 1
fi

# ============================================================================ #
# Generate mesh and run case

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 4 0 1 0 1 $Nx $Ny $Nz .false. .false. .false.

for case in $cases; do
    echo "Running case: $case"
    # neko $case >${case%.*}.log
done

# End of file
# ============================================================================ #
