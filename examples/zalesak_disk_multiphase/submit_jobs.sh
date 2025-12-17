#!/bin/bash

# ============================================
# Job submission script for parameter sweep
# ============================================

# ============================================
# CHOOSE MODE: Uncomment one of the following
# ============================================

# MODE 1: Test with single job
TEST_DIR="gamma_0.01_epsilon_0.02"  # Specify which combination to test
MODE="single"

# MODE 2: Submit all as job array (uncomment for full sweep)
# MODE="array"

# ============================================

if [ "$MODE" == "single" ]; then
    # Submit single test job
    echo "Submitting single test job for: $TEST_DIR"

    sbatch <<EOF
#!/bin/bash
#SBATCH -A naiss2025-22-129
#SBATCH -t 03:00:00
#SBATCH -N 1
#SBATCH -J test_${TEST_DIR}
#SBATCH -p main

# Load environment modules
module unload PDC
ml PrgEnv-cray/8.5.0
ml cpe/23.09

# Source environment file if it exists
if [ -f \$SLURM_SUBMIT_DIR/.env ]; then
    source \$SLURM_SUBMIT_DIR/.env
fi

cd \$SLURM_SUBMIT_DIR/${TEST_DIR}
srun -u -n 64 ../neko zalesak.case
EOF

elif [ "$MODE" == "array" ]; then
    # Submit job array for all combinations

    # Get list of all directories
    DIRS=(gamma_*)
    N=${#DIRS[@]}

    echo "Submitting job array with $N jobs..."

    # Write directory list to file
    printf '%s\n' "${DIRS[@]}" > sweep_dirs.txt

    sbatch <<EOF
#!/bin/bash
#SBATCH -A naiss2025-22-129
#SBATCH -t 03:00:00
#SBATCH -N 1
#SBATCH -J sweep_array
#SBATCH -p main
#SBATCH --array=0-$((N-1))

# Load environment modules
module unload PDC
ml PrgEnv-cray/8.5.0
ml cpe/23.09

# Source environment file if it exists
if [ -f \$SLURM_SUBMIT_DIR/.env ]; then
    source \$SLURM_SUBMIT_DIR/.env
fi

# Get directory for this array task
DIR=\$(sed -n "\$((SLURM_ARRAY_TASK_ID + 1))p" sweep_dirs.txt)

echo "Running job \$SLURM_ARRAY_TASK_ID: \$DIR"

cd \$SLURM_SUBMIT_DIR/\$DIR
srun -u -n 64 ../neko zalesak.case
EOF

else
    echo "Error: MODE must be 'single' or 'array'"
    exit 1
fi
