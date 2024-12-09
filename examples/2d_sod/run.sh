# Generate mesh
genmeshbox 0 1 0 1 0 0.1 30 30 1 .true .true. .true.

# Run
mpirun -n 12 neko sod.case