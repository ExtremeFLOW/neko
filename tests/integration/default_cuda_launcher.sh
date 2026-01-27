#!/bin/bash

# Usage: ./run_neko.sh <nprocs> <case_file> <neko_path>

# Get the input parameters
nprocs="$1"
case_file="$2"
neko="$3"

# Check that all three parameters are provided
if [ -z "$nprocs" ] || [ -z "$case_file" ] || [ -z "$neko" ]; then
    echo "Usage: $0 <nprocs> <case_file> <neko_path>"
    exit 1
fi

# Run the command
mpirun -np "$nprocs" --bind-to core --map-by ppr:1:core \
    bash -c "CUDA_VISIBLE_DEVICES=\$((OMPI_COMM_WORLD_RANK)) $neko $case_file"

