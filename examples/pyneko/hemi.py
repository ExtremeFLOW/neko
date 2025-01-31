import pyneko
import json

# Initialise Neko
pyneko.init()
pyneko.job_info()

# Create a Neko case from a JSON file
hemi_json = json.load(open('hemi.case'))
hemi_case = pyneko.case_init(hemi_json)

# To solve the entire case we can call solve()
# pyneko.solve(hemi_case

# To manually step forward in time, call step()
t = 0.0
for tstep in range (1, 3):
    pyneko.step(hemi_case, t, tstep)
    t = t + hemi_json['case']['timestep']
    pyneko.output(hemi_case, t, tstep)

# Force an output
pyneko.output(hemi_case, t, tstep, force=True)

# Cleanup
pyneko.case_free(hemi_case)
pyneko.finalize()
