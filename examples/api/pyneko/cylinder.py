import pyneko
import json

# Initialise Neko
pyneko.init()
pyneko.job_info()

# Create a Neko case from a JSON file
cylinder_json = json.load(open('cylinder.case'))
cylinder_case = pyneko.case_init(cylinder_json)

# To solve the entire case we can call solve()
#pyneko.solve(cylinder_case)

# To manually step forward in time, call step()
while pyneko.time(cylinder_case) < pyneko.end_time(cylinder_case):
    pyneko.step(cylinder_case)

# To retrive a dataclass representation of a field with a numpy array
# of the data in a field from Neko's field registry use field()
u = pyneko.field(b'u')

# Cleanup
pyneko.case_free(cylinder_case)
pyneko.finalize()
