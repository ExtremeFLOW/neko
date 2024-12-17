# Generate mesh
genmeshbox 0 1 0 0.5 0 0.01 100 50 1 .false. .true. .true.
# genmeshbox 0 1 0 0.01 0 0.001 1000 10 1 .false. .true. .true.
# genmeshbox 0 1 0 0.001 0 0.0001 10000 10 1 .false. .true. .true.

# Run
mpirun -n 12 neko sod.case