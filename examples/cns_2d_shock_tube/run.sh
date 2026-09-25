#!/bin/bash
# Generate 2D mesh: [0,1] x [0,1], 50x50 elements, periodic only in z
genmeshbox 0 1 0 1 0 0.02 50 50 1 .false. .false. .true.

# Compile user file
makeneko 2d_shock_tube.f90

# Run
mpirun -np 4 ./neko 2d_shock_tube.case
