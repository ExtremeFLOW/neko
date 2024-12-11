# Generate mesh
genmeshbox 0 1 0 1 0 0.03 30 30 1 .false. .true. .true.

# Run
mpirun -n 12 neko sod.case