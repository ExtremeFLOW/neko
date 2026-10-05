"""Emit one .case per (H, N, nsvv_ratio) cell of Saini section 4.4.

Usage: gen_cases.py <outdir> <eps-rule> <Hden>:<N>[:<ratio>] ...

<eps-rule> is mandatory and has no default. It is the width of the sign function,
Eq. (46). `0.25` (absolute) is what the committed cells carry: the authors' value,
from their code (`signls`; README.md section 2). `HN:1` is the reading of the
paper's Eq. (36), eps = xi*H with their xi = 1/N. Forms:

    HN:1         C * H/N        (the paper's Eq. 36 at C = 1)
    H:4.0        C * H
    hmax:8.9     C * h_GLL,max  = C * HMAX[N] * H
    0.4          absolute length

Cells differ only in mesh_file, polynomial_order, epsilon and the SVV ratio;
everything else is fixed by the paper.
"""
import json, os, sys

MESH = {5: "box20.nmsh", 10: "box40.nmsh", 20: "box80.nmsh"}
# GLL spacing in units of H, max over the element (the centre gap).
HMAX = {3: 0.4472, 4: 0.3273, 5: 0.2852, 6: 0.2344, 7: 0.2093, 8: 0.1816}

def epsilon(rule, H, n):
    if ":" not in rule:
        return float(rule)
    kind, c = rule.split(":", 1)
    c = float(c)
    return {"hmax": c*HMAX[n]*H, "H": c*H, "HN": c*H/n}[kind]

def case(hden, n, ratio, eps_rule):
    H = 1.0 / hden
    return {
        "version": 1.0,
        "case": {
            "mesh_file": MESH[hden],
            "output_directory": "output_%s" % name(hden, n, ratio),
            "output_boundary": False,
            "output_checkpoints": False,
            "output_at_end": False,
            "output_precision": "double",
            "time": {
                "end_time": 0.0,
                "variable_timestep": False,
                "timestep": 1e-06
            },
            "numerics": {
                "time_order": 3,
                "polynomial_order": n,
                "dealias": True
            },
            "fluid": {
                "scheme": "pnpn",
                "freeze": True,
                "initial_condition": {"type": "uniform", "value": [0.0, 0.0, 0.0]},
                "velocity_solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "projection_space_size": 0,
                    "absolute_tolerance": 1e-06,
                    "max_iterations": 800
                },
                "pressure_solver": {
                    "type": "gmres",
                    "preconditioner": {"type": "hsmg"},
                    "projection_space_size": 0,
                    "absolute_tolerance": 0.001,
                    "max_iterations": 800
                },
                "output_control": "simulationtime",
                "output_value": 1.0,
                "boundary_conditions": [
                    {"type": "velocity_value",
                     "zone_indices": [1, 2, 3, 4],
                     "value": [0.0, 0.0, 0.0]}
                ]
            },
            "cdi": {
                "epsilon": epsilon(eps_rule, H, n),
                "report_every": 500,
                "redistance": {
                    "dtau": 0.0005,
                    "tau_end": 6.0,
                    "svv": {"c0": 2.0, "nsvv_ratio": float(ratio)}
                }
            },
            "scalars": [
                {
                    "name": "psi",
                    "enabled": True,
                    "cp": 1.0,
                    "lambda": 1e-16,
                    "initial_condition": {"type": "user"},
                    "solver": {
                        "type": "cg",
                        "preconditioner": {"type": "jacobi"},
                        "projection_space_size": 0,
                        "absolute_tolerance": 1e-06,
                        "max_iterations": 800
                    },
                    "boundary_conditions": [
                        {"type": "neumann",
                         "zone_indices": [1, 2, 3, 4],
                         "flux": 0.0}
                    ]
                }
            ]
        }
    }

def name(hden, n, ratio):
    base = "circles_h%d_n%d" % (hden, n)
    return base if ratio == 6 else base + "_nsvv%d" % ratio

outdir, eps_rule = sys.argv[1], sys.argv[2]
for spec in sys.argv[3:]:
    p = spec.split(":")
    hden, n = int(p[0]), int(p[1])
    ratio = int(p[2]) if len(p) > 2 else 6
    fn = os.path.join(outdir, name(hden, n, ratio) + ".case")
    with open(fn, "w") as f:
        json.dump(case(hden, n, ratio, eps_rule), f, indent=4)
        f.write("\n")
    print(fn, "eps = %.6g" % epsilon(eps_rule, 1.0/hden, n))
