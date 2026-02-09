#!/bin/bash

# Compile the Zalesak case
makeneko zalesak.f90
chmod +x neko

# Array of gamma values to test
gammas=("0" "0.000001" "0.00001" "0.0001" "0.001" "0.002" "0.005" "0.006" "0.007" "0.008" "0.009" "0.01" "0.02" "0.03" "0.04" "0.05" "0.1" "0.2" "0.5" "1.0")

# Array of epsilon values to test
epsilon=("0.005" "0.01" "0.015" "0.02" "0.025" "0.03" "0.035" "0.04" "0.045" "0.05")

# make a directory for the specific epsilon gamma runs
mkdir -p zalesak_sweep
cd zalesak_sweep



# Loop over each gamma value
for g in "${gammas[@]}"; do
    for e in "${epsilon[@]}"; do
        dir="gamma_${g}_epsilon_${e}"
        mkdir -p "$dir"
        cd "$dir"

        # Make a template case file
        sed "s/GAMMA_VALUE/$g/g" template.case > "gamma_$g.case"

        # make sure the mesh can be reached from this directory

        # Return to the sweep directory
        cd ..
    done

    # Create a directory for the current gamma-epsilon value
