import json

def create_case(mesh_file="", lx=4, nsamples=0, T_end=0.0):
    default = {
        "mesh_file" : mesh_file,
        "fluid_scheme" : "pnpn",
        "lx" : lx,
        "nsamples" :  nsamples,
        "dt" : 0.001,
        "T_end" : T_end,
        "ksp_vel" : {
            "type" : "cg",
            "pc" : "jacobi",
            "abstol" : 1e-09
            
        },
        "ksp_prs": {
            "type" : "gmres",
            "pc" : "hsmg",
            "abstol" : 1e-09
        }
    }
    return json.loads(json.dumps(default))


