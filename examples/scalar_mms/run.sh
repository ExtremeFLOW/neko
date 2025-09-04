# Compile the case
makeneko scalar_mms.f90

# Generate mesh
genmeshbox 0 1 0 1 0 1 20 1 1 .false. .true. .true.
# Run the code
./neko scalar_mms.case
