#!/usr/bin/bash

# ============================================================================ #
# Define the help function

function help() {
    echo -e "prepare.sh case"
    echo -e "  Generate a mesh to prepare for the desired case."
    echo -e "  The input arguments are the number of cells in the x, y, and z"
    echo -e "  directions, respectively."
    echo -e ""
    echo -e "  If no input arguments are provided, the default mesh size is"
    echo -e "  32x8x8."
    echo -e ""
    echo -e "  Example usage:"
    echo -e "    prepare.sh -x32 -y8 -z8"
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    echo -e "  -x#         Number of cells in the x direction."
    echo -e "  -y#         Number of cells in the y direction."
    echo -e "  -z#         Number of cells in the z direction."
    echo -e "  -q, --quiet Suppress output."
    echo -e ""
    echo -e "  See Readme for additional details."
    exit 0
}

# Handle options
Nx=32 && Ny=8 && Nz=8
for arg in "$@"; do
    if [ "${arg:0:2}" == "--" ]; then
        case ${arg:2} in
        help) help ;;
        quiet) QUIET=1 ;;
        *) echo -e "Invalid option: $arg" >&2 && help ;;
        esac
    elif [ "${arg:0:1}" == "-" ]; then
        case ${arg:1:1} in
        h) help ;;
        x) Nx=${arg:2} ;;
        y) Ny=${arg:2} ;;
        z) Nz=${arg:2} ;;
        q) QUIET=1 ;;
        *) echo -e "Invalid option: ${arg:1}" >&2 && help ;;
        esac
    else
        cases="$cases $arg"
    fi
done

# ============================================================================ #
# Ensure Neko can be found and set default mesh size

if [ "$NEKO_DIR" ]; then
    PATH=$NEKO_DIR/bin:$PATH
fi

if [[ -z $(which genmeshbox) ]]; then
    echo -e "Neko not found." >&2
    echo -e "Please ensure Neko is installed and in your PATH." >&2
    echo -e "Alternatively, set the NEKO_DIR environment variable." >&2
    exit 1
fi

# ============================================================================ #
# Generate mesh and prepare case

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 4 0 1 0 1 $Nx $Ny $Nz .false. .false. .false.

# End of file
# ============================================================================ #
