# Convert mesh
rea2nbin hemi.re2

# Run
mpirun -n 4 neko hemi.case
