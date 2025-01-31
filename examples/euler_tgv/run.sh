#!/bin/bash
# Generate mesh
genmeshbox 0 6.28318530718 0 6.28318530718 0 6.28318530718 20 20 20 .true. .true. .true.

# Run
mpirun -n 4 neko euler_tgv.case