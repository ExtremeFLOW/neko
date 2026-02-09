#!/bin/bash

# ============================================
# Job submission script for postprocessing
# Advecting drop multiphase test case
# ============================================

# Directory containing simulation data
DATA_DIR="/cfs/klemming/scratch/e/eriksie/advecting_drop_param_analysis_data"

# Output directory for plots
OUTPUT_DIR="/cfs/klemming/projects/supr/kthmech/eriksie/data/advecting_drop_reference_plots"

sbatch <<EOF
#!/bin/bash
#SBATCH -A naiss2025-22-129
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J postprocess_advdrop
#SBATCH -p shared

# Load Python module
module load cray-python

# Activate pysemtools virtual environment
source /cfs/klemming/projects/supr/kthmech/eriksie/pysemtools_venv/bin/activate

cd \$SLURM_SUBMIT_DIR

echo "Running reference-style analysis..."
echo "Data directory: ${DATA_DIR}"
echo "Output directory: ${OUTPUT_DIR}"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Run the reference-style analysis (--base-dir processes all gamma_* subdirs)
mpirun -n 1 python analyze_reference_style.py --base-dir ${DATA_DIR} --output-dir ${OUTPUT_DIR}

echo "Postprocessing complete!"
EOF

echo "Submitted postprocessing job for: ${DATA_DIR}"
