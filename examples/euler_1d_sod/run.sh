#!/bin/bash
# Generate mesh
genmeshbox 0 1 0 0.01 0 0.01 100 1 1 .false. .true. .true.
# genmeshbox 0 1 0 0.001 0 0.001 1000 1 1 .false. .true. .true.
# genmeshbox 0 1 0 0.0001 0 0.0001 10000 1 1 .false. .true. .true.

# Run
mpirun -n 4 ./neko sod.case
