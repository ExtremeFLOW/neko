#!/bin/bash
# linting helper script for Fortran files
# This script uses the flinter tool to lint the Fortran files
# The linter tool can be installed using the following command:
# pip install flinter nobvisual

# ---------------------------------------------------------------------------- #
# Helper function

function help() {
    echo "Usage: lint.sh [options] [files]"
    echo "Lint Fortran files in the repository."
    echo "This tool will lint the Fortran files using the flinter tool."
    echo "It can be used for specific files or all the modified files in the "
    echo "repository. Relative paths will be relative to the Neko root directory."
    echo ""
    echo "Options:"
    echo "  -h, --help                Print this help message"
    echo "      --target_branch       The branch to compare with (default: develop)"
    echo ""
    echo "Arguments:"
    echo "  files                     List of files to lint"
    echo ""
    echo "If no files are provided, lint all the Fortran files in the repository which are modified."
}

# ---------------------------------------------------------------------------- #
# Handle options and arguments

# Set the root directory
ROOT_DIR=$(realpath $(dirname $0)/../../)

# Check if the flint tool is installed
if ! command -v flint &>/dev/null; then
    echo "The flint tool is not installed. Please install it using the following command:"
    echo "pip install flint nobvisual"
    exit 1
fi

# List possible options
OPTIONS=help,target_branch:
OPT=h

# Default values for the options
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
    "-h" | "--help") help && exit ;;                    # Print help
    "--target_branch") TARGET_BRANCH="$2" && shift 2 ;; # Target branch

    # End of options
    "--") shift && break ;;
    *) echo "Unknown option: $1" >&2 && help && exit 1 ;; # Unknown option
    esac
done

# ---------------------------------------------------------------------------- #
# Handle the files

# If no files are provided, lint all the Fortran files in the repository which
# are modified.
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

# ---------------------------------------------------------------------------- #
# Lint the files

failed_files=()
printf "Linting files:\n"
for file in ${files[@]}; do

    # If the file is not a Fortran file, skip it.
    if [[ ${file: -4} != ".f90" && ${file: -4} != ".F90" ]]; then
        continue
    fi

    printf "\t- $file"
    score=$(flint score -r $ROOT_DIR/flinter_rc.yml $file 2>/dev/null |
        sed -n 's/.*|\([^|]*\)|.*/\1/p')
    printf ": $score\n"

    if (($(echo "$score < 10" | bc -l))); then
        failed_files+=($file)
    fi
done

if [ ${#failed_files[@]} -eq 0 ]; then
    echo "All files are linted successfully."
    exit 0
fi

# ---------------------------------------------------------------------------- #
# Print possible improvements

printf "Files that can be improved:\n" | tee linter-report.txt
for file in ${failed_files[@]}; do
    printf "\t- $file\n" | tee -a linter-report.txt
done

if [ ${#failed_files[@]} -gt 0 ]; then
    for file in ${failed_files[@]}; do
        report=$(flint lint -r $ROOT_DIR/flinter_rc.yml $file)
        if [ -z "$report" ]; then
            report=$(flint stats -r $ROOT_DIR/flinter_rc.yml $file)
        fi

        printf "%.s-" {1..80} | tee -a linter-report.txt
        printf "\n" | tee -a linter-report.txt
        printf "Linting improvements for \n\t$file\n\n" | tee -a linter-report.txt
        echo "$report" | tee -a linter-report.txt
    done
fi
