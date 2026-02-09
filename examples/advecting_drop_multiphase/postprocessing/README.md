# Postprocessing Scripts

Quick reference for analyzing advecting droplet simulation results.

**Prerequisites**: Python 3.x with `pysemtools`, `numpy`, `matplotlib`, `mpi4py`

---

## Quick Start

```bash
cd postprocessing

# Check initial/final conditions
python3 check_initial_condition.py ../case_directory/field0.f00000 <epsilon>

# Create animation
python3 create_animation.py --mode animation \
    --output-dir ../case_directory \
    --animation-file ../case_directory/animation.mp4

# Analyze parameter sweep
python3 analyze_reference_style.py
```

---

## Scripts

### check_initial_condition.py

Verify field matches analytical 2D tanh profile

```bash
python3 check_initial_condition.py <field_file> <epsilon>
```

### create_animation.py

Create MP4 video or PNG frames from field files

```bash
# MP4 video
python3 create_animation.py --mode animation \
    --output-dir ../gamma_0.01_epsilon_0.02 \
    --animation-file advecting_drop.mp4

# Individual frames
python3 create_animation.py --mode frames \
    --output-dir ../gamma_0.01_epsilon_0.02 \
    --frames-dir animation_frames
```

### analyze_reference_style.py

Generate heatmaps and line plots for parameter sweep (saves to `reference_style_plots/`)

```bash
python3 analyze_reference_style.py
```

### debug_slice.py

Debug 2D slice extraction at z≈0

```bash
python3 debug_slice.py
```
