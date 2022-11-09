import json

def create_case(mesh_file="", 
                lx=4, 
                source_term="noforce",
                initial_condition="uniform",
                nsamples=0, 
                dt=0.001,                
                T_end=0.0,
                uinf=[1.0,0.0,0.0]) :
    default = {
        "case" : {
            "mesh_file" : mesh_file,
            "fluid_scheme" : "pnpn",
            "lx" : lx,
            "source_term" : source_term,
            "initial_condition" : initial_condition
        },
        "parameters" : {
            "nsamples" :  nsamples,        
            "dt" : dt,
            "T_end" : T_end,
            "uinf" : uinf,
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
    }
    return json.loads(json.dumps(default))


