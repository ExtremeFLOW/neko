# Compile the case
makeneko advecting_cone.f90
# Generate probes
python create_probes.py

# Remove results if exist
[ -f output.csv ] && rm output.csv

# Set order
python change_order.py 2
# Generate mesh
genmeshbox -2 2 -2 2 0 0.1 30 30 1 .true .true. .true.
# Run the code
mpirun -n 8 ./neko advecting_cone.case
# Save result
mv output.csv output_2.csv

# Order 3
python change_order.py 3
genmeshbox -2 2 -2 2 0 0.1 22 22 1 .true .true. .true.
mpirun -n 8 ./neko advecting_cone.case
mv output.csv output_3.csv

# Order 5
python change_order.py 5
genmeshbox -2 2 -2 2 0 0.1 15 15 1 .true .true. .true.
mpirun -n 8 ./neko advecting_cone.case
mv output.csv output_5.csv
