# Bound Violation Analysis

This directory contains scripts for post-processing bound violations in the Zalesak disk parameter sweep.

## Files

- **`analyze_bound_violations.py`**: Main script for sweeping through all cases and analyzing bound violations
- **`postprocess.ipynb`**: Original Jupyter notebook with various analysis methods

## Requirements

The analysis script requires:
- Python 3.x
- numpy
- matplotlib
- mpi4py
- pysemtools (for reading Neko field files)

These should be the same dependencies as the Jupyter notebook.

## Usage

### Basic Usage

From the `postprocessing` directory, run:

```bash
python3 analyze_bound_violations.py
```

Or using MPI:

```bash
mpirun -n 1 python3 analyze_bound_violations.py
```

### What the Script Does

The script will:

1. **Find all cases**: Searches for directories matching `gamma_*_epsilon_*` in the parent directory
2. **Read field data**: Loads the final time step (`field0.f00002`) from each case
3. **Compute violations**: Calculates bound violations for the scalar field (expected to be in [0, 1])
4. **Generate plots**:
   - Maximum violation vs γ (for each ε)
   - Maximum violation vs ε (for each γ)
   - 2D heatmap of violations (γ-ε parameter space)
   - Violation fraction vs γ (percentage of points violating bounds)
5. **Print statistics**: Shows worst/best cases and summary information

### Output Files

The script generates the following plots in the `postprocessing` directory:

- `bound_violation_vs_gamma.png` - Line plots showing max violation vs γ for different ε values
- `bound_violation_vs_epsilon.png` - Line plots showing max violation vs ε for different γ values
- `bound_violation_heatmap.png` - 2D heatmap showing violations across γ-ε parameter space
- `violation_fraction_vs_gamma.png` - Percentage of points violating bounds vs γ

## Understanding Bound Violations

The scalar field (temperature/phase field) should theoretically remain bounded in [0, 1]. Violations occur when:
- **Upper violations**: Field value > 1.0
- **Lower violations**: Field value < 0.0

The script computes:
- **Maximum violation**: Largest absolute violation magnitude across all points
- **Mean violation**: Average violation magnitude
- **Violation count**: Number of points violating bounds
- **Violation fraction**: Percentage of total points with violations

## Example Output

```
================================================================================
Bound Violation Analysis for Zalesak Disk Parameter Sweep
================================================================================

Found 200 case directories

✓ gamma_0.01_epsilon_0.02              | γ=0.01     ε=0.02     | Max violation: 0.000123 | Violations: 45/32768
...

================================================================================
Successfully processed 200 cases
================================================================================

Unique gamma values: 20
Unique epsilon values: 10

Gamma range: [0.0, 1.0]
Epsilon range: [0.005, 0.05]

Max violation range: [0.000000, 0.123456]
```

## Customization

You can modify the script to:
- Change the field file to analyze (default: `field0.f00002`)
- Adjust bound limits (default: [0, 1])
- Add additional statistics or plots
- Change output directory

Edit the configuration section in `main()`:

```python
# Configuration
base_dir = Path('..')  # Parent directory containing gamma_* folders
output_dir = Path('.')  # Current directory for outputs
field_file_name = 'field0.f00002'  # Final time step
```

## Troubleshooting

**Q: "No case directories found"**
A: Make sure you're running from the `postprocessing` directory and that `gamma_*_epsilon_*` directories exist in the parent directory.

**Q: "field0.f00002 not found"**
A: The simulation may not have completed. Check that the field files exist in each case directory.

**Q: Import errors for pysemtools**
A: Ensure pysemtools is installed and you're using the correct Python environment.

## Next Steps

After running the analysis:
1. Review the generated plots to identify parameter regions with low violations
2. Compare violation patterns with other metrics (shape error, boundary error, etc.)
3. Use the summary statistics to select optimal γ and ε values
4. Consider running additional cases in parameter regions of interest
