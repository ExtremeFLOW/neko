import pyNeko

# Initialize Neko
pyNeko.init()

# Create a json case
hemi_case = pyNeko.create_case(mesh_file="hemi.nmsh", lx=6, T_end=1e-3)

# Solve the case
pyNeko.solve(hemi_case)

# Increase sample rate and rerun to get some data
hemi_case['parameters']['nsamples'] = 1
pyNeko.solve(hemi_case)
