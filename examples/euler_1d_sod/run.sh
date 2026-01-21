#!/bin/bash
# Generate mesh with domain [0,1]
# Discontinuity at x=0.5
genmeshbox 0 1 0 0.001 0 0.001 1000 1 1 .false. .true. .true.

# Run
mpirun -n 4 ./neko sod.case
