# Generate mesh
genmeshbox 0 1 0 1 0 0.01 100 100 1 .true. .true. .true.

# Run
mpirun -n 12 neko compressible_smooth.case
