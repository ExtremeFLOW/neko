#!/bin/bash

# ============================================
# Setup script for gamma-epsilon parameter sweep
# ============================================

# Paths for Dardel
MAKENEKO="/cfs/klemming/projects/supr/naiss2025-22-129/eriksie/software/neko-cpu/bin/makeneko"
OUTPUT_BASE="/cfs/klemming/scratch/e/eriksie/zalesak_param_analysis_data"

# ============================================
# MODIFY THESE ARRAYS FOR TESTING vs FULL SWEEP
# ============================================

# For testing - use single values:
gammas=("0.01")
epsilons=("0.02")

# For full sweep - uncomment and use:
# gammas=("0" "0.000001" "0.00001" "0.0001" "0.001" "0.002" "0.005" "0.006" "0.007" "0.008" "0.009" "0.01" "0.02" "0.03" "0.04" "0.05" "0.1" "0.2" "0.5" "1.0")
# epsilons=("0.005" "0.01" "0.015" "0.02" "0.025" "0.03" "0.035" "0.04" "0.045" "0.05")

# ============================================

echo "Compiling Neko with zalesak.f90..."
$MAKENEKO zalesak.f90
chmod +x neko

echo "Creating sweep directories..."

for gamma in "${gammas[@]}"; do
    for epsilon in "${epsilons[@]}"; do
        DIR="gamma_${gamma}_epsilon_${epsilon}"
        echo "Setting up: $DIR"

        # Create directory
        mkdir -p "$DIR"

        # Create output directory on scratch
        OUTPUT_DIR="${OUTPUT_BASE}/gamma_${gamma}_epsilon_${epsilon}"

        # Copy template and modify epsilon, gamma, output_directory, and mesh path
        sed -e "s/\"epsilon\": 0.01/\"epsilon\": ${epsilon}/" \
            -e "s/\"gamma\": 0.01/\"gamma\": ${gamma}/" \
            -e "s|\"output_directory\":\".*\"|\"output_directory\":\"${OUTPUT_DIR}\"|" \
            -e "s|\"mesh_file\": \"40x40.nmsh\"|\"mesh_file\": \"../40x40.nmsh\"|" \
            zalesak_sweep.case > "${DIR}/zalesak.case"
    done
done

echo "Setup complete!"
echo "Total combinations: $((${#gammas[@]} * ${#epsilons[@]}))"
