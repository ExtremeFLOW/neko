#!/bin/bash

# ============================================
# Setup script for gamma-epsilon parameter sweep
# ============================================

# Paths for Dardel
MAKENEKO="/cfs/klemming/projects/supr/naiss2025-22-129/eriksie/software/neko-cpu/bin/makeneko"
OUTPUT_BASE="/cfs/klemming/scratch/e/eriksie/advecting_drop_param_analysis_data"

# ============================================
# MODIFY THESE ARRAYS FOR TESTING vs FULL SWEEP
# ============================================

# For testing - use single values:
gammas=("0.01")
epsilons=("0.02")

# For full sweep - uncomment and use:
# Updated with additional values to refine parameter space
# Gamma: Added 0.003, 0.004, 0.011, 0.012, 0.013, 0.015, 0.06, 0.07, 0.15 (29 values total)
# gammas=("0" "0.000001" "0.00001" "0.0001" "0.001" "0.002" "0.003" "0.004" "0.005" "0.006" "0.007" "0.008" "0.009" "0.01" "0.011" "0.012" "0.013" "0.015" "0.02" "0.03" "0.04" "0.05" "0.06" "0.07" "0.1" "0.15" "0.2" "0.5" "1.0")
# Epsilon: Added 0.0025, 0.055, 0.06 (13 values total)
# epsilons=("0.0025" "0.005" "0.01" "0.015" "0.02" "0.025" "0.03" "0.035" "0.04" "0.045" "0.05" "0.055" "0.06")

# ============================================

echo "Compiling Neko with advecting_drop.f90..."
$MAKENEKO advecting_drop.f90
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
            -e "s|\"mesh_file\": \"box.nmsh\"|\"mesh_file\": \"../box.nmsh\"|" \
            advecting_drop_sweep.case > "${DIR}/advecting_drop.case"
    done
done

echo "Setup complete!"
echo "Total combinations: $((${#gammas[@]} * ${#epsilons[@]}))"
