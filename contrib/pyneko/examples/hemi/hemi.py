import pyNeko
import json

# Initialize Neko
pyNeko.init()

# Create a json case
hemi_case = json.load(open('hemi.case'))

# Solve the case
pyNeko.solve(hemi_case)

# Increase sample rate and rerun to get more data
hemi_case['case']['nsamples'] = 2
pyNeko.solve(hemi_case)


# Finalize Neko
pyNeko.finalize()
