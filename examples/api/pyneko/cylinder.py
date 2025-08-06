import pyneko
import json
import ctypes

# Initialise Neko
pyneko.init()
pyneko.job_info()


cylinder_json = json.load(open('cylinder.case'))

# Define initial conditions
def initial(name, len):
    scheme_name = ctypes.string_at(name, len).decode()

    # In this case the check is redudnant, but if both fluid and scalar
    # initial conditions should be defined, `scheme_name` can be used
    # to identifiy which solver has called the callback
    if (scheme_name == "fluid"):
        # To retrive a dataclass representation of a field with a numpy array
        # of the data in a field from Neko's field registry use field()
        u = pyneko.field(b'u')
        for i in range(0, u.dm.ntot):
            u.x[i] = 1.0

# Define inflow conditions
def inflow(msk, msk_size, t, tstep):
    bc_mask = pyneko.bc_mask(msk, msk_size)

    # In this case this check is redudant, but if different user provied
    # boundary conditions share the same callback, one could use the content
    # of the callbacks field list to identify which condition should be applied,
    # for example a velocity condition passes the fields (u, v, w)
    if (pyneko.callback_field_name(1, b'u')): # Note: Fortran indices
        u = pyneko.callback_field(b'u')
        v = pyneko.callback_field(b'v')
        w = pyneko.callback_field(b'w')
        for i in range(0, msk_size):
            idx = bc_mask[i] - 1 # Neko's mask is based on Fortran indices
            u.x[idx] = 1.0
            v.x[idx] = 0.0
            w.x[idx] = 0.0

# Create a Neko callback for the initial and inflow conditions
cb_cylinder_ic = pyneko.initial_condition(initial)
cb_cylinder_if = pyneko.dirichlet_condition(inflow)

# Create a Neko case from a JSON file and provied (optional) callbacks
cylinder_case = pyneko.case_init(cylinder_json,
                                 cb_initial_condition = cb_cylinder_ic,
                                 cb_dirichlet_condition = cb_cylinder_if)

# To solve the entire case we can call solve()
#pyneko.solve(cylinder_case)

# To manually step forward in time, call step()
while pyneko.time(cylinder_case) < pyneko.end_time(cylinder_case):
    pyneko.step(cylinder_case)

# Cleanup
pyneko.case_free(cylinder_case)
pyneko.finalize()
