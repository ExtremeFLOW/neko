#!/bin/bash
# Formatting helper script for Fortran files
# This script uses the findent tool to format the Fortran files
# The findent tool can be installed using the following command:
# pip install findent

# ---------------------------------------------------------------------------- #
# Helper function

function help() {
    echo "Usage: format.sh [options] [files]"
    echo "Format Fortran files in the repository."
    echo "This tool will format the Fortran files using the findent tool."
    echo "It can be used for specific files or all the modified files in the "
    echo "repository. Relative paths will be relative to the Neko root directory."
    echo ""
    echo "Options:"
    echo "  -h, --help                Print this help message"
    echo "  -k, --indent_continuation Set the continuation indentation (default: 5)"
    echo "      --target_branch       The branch to compare with (default: develop)"
    echo ""
    echo "Arguments:"
    echo "  files                     List of files to format"
    echo ""
    echo "If no files are provided, format all the Fortran files in the repository which are modified."
}

# ---------------------------------------------------------------------------- #
# Handle options and arguments

# Set the root directory
ROOT_DIR=$(realpath $(dirname $0)/../../)

# Check if the findent tool is installed
if ! command -v findent &>/dev/null; then
    echo "The findent tool is not installed. Please install it using the following command:"
    echo "pip install findent"
    exit 1
fi

# List possible options
OPTIONS=help,indent_continuation:,target_branch:
OPT=h,k:

# Assign default values to the options
CONTINUATION="5"
TARGET_BRANCH="develop"

# Parse the inputs for options
if [ $(uname) == "Darwin" ]; then
    # Use gnu-getopt on macOS
    if command -v gnu-getopt >/dev/null 2>&1; then
        PARSED=$(gnu-getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
    else
        echo "GNU-getopt not found." >&2
        echo "Falling back to BSD getopt, long options are not supported." >&2
        PARSED=$(getopt $OPT "$@")
    fi
else
    # Use getopt on Linux
    PARSED=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
fi
eval set -- "$PARSED"

# Loop through the options and set the variables
while true; do
    case "$1" in
    "-h" | "--help") help && exit ;;                                # Print help
    "-k" | "--indent_continuation") CONTINUATION="$2" && shift 2 ;; # Continuation indentation
    "--target_branch") TARGET_BRANCH="$2" && shift 2 ;;             # Target branch

    # End of options
    "--") shift && break ;;
    *) echo "Unknown option: $1" >&2 && help && exit 1 ;; # Unknown option
    esac
done

[ $CONTINUATION == "none" ] && CONTINUATION="-"

# Set the options
FLAGS="-Rr -i2 -d3 -f3 -s3 -c3 -w3 -t3 -j3 --ws_remred --openmp=0"
FLAGS="$FLAGS -k$CONTINUATION"

# ---------------------------------------------------------------------------- #
# Handle the files

# If no files are provided, format all the Fortran files in the repository which are modified
if [ $# -gt 0 ]; then
    files=($@)
else
    echo "Formatting all the modified Fortran files in the repository."

    git fetch origin $TARGET_BRANCH >/dev/null 2>&1
    files=($(git diff --name-only --diff-filter=d origin/$TARGET_BRANCH))
fi

# Prepend the root directory to the files if they are not absolute paths and
# trim the leading ./ if present
for i in ${!files[@]}; do
    if [[ ${files[$i]} != /* ]]; then
        files[$i]="$(realpath $ROOT_DIR/${files[$i]})"
    fi
done

if [ ${#files[@]} -eq 0 ]; then
    echo "No files to format"
    exit 0
fi

# ---------------------------------------------------------------------------- #
# Format the files

for file in ${files[@]}; do
    if [[ ${file: -4} != ".f90" && ${file: -4} != ".F90" ]]; then
        continue
    fi

    printf "\t- $file\n"
    findent $FLAGS <$file >$file.tmp
    mv -f $file.tmp $file
done
