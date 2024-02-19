# Compile the case
makeneko scalar_mms.f90
# Generate probes
python create_probes.py

# Remove results if exist
rm output.csv

# Generate mesh
genmeshbox 0 1 0 1 0 0.1 10 1 1 .true .true. .true.
# Run the code
./neko scalar_mms.case
