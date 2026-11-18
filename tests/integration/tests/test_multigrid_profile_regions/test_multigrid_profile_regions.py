"""Regression test for runtime statistics with both multigrid preconditioners.

"""
from os.path import join
import json

import json5
from testlib import get_neko, run_neko, configure_nprocs, get_neko_dir


def test_multigrid_profile_regions(launcher_script, log_file, tmp_path):
    """The hsmg and phmg preconditioners used in one case each get their own
    profiling region in the runtime statistics."""
    neko_dir = get_neko_dir()
    with open(join(neko_dir, "examples", "cylinder", "cylinder.case")) as f:
        case = json5.load(f)

    case_object = case["case"]
    case_object["mesh_file"] = join("meshes", "small_test_cyl.nmsh")
    case_object["output_directory"] = str(tmp_path)
    case_object["output_at_end"] = False
    case_object["output_boundary"] = False
    case_object["output_checkpoints"] = False
    case_object["time"] = {"end_time": 2.0e-3, "timestep": 1.0e-3}
    case_object["runtime_statistics"] = {"enabled": True}
    case_object.pop("simulation_components")

    fluid = case_object["fluid"]
    fluid["output_control"] = "never"
    # The pressure keeps the hsmg preconditioner of the example
    assert fluid["pressure_solver"]["preconditioner"]["type"] == "hsmg"

    # A passive scalar takes phmg. Only applying it matters here, not
    # converging with it, so its iterations are capped to keep the test short.
    case_object["scalar"] = {
        "enabled": True,
        "lambda": 0.007,
        "cp": 1.0,
        "initial_condition": {"type": "uniform", "value": 0.0},
        "boundary_conditions": [
            {"type": "dirichlet", "value": 1.0, "zone_indices": [7]},
            {"type": "dirichlet", "value": 0.0, "zone_indices": [1]},
            {"type": "neumann", "flux": 0.0, "zone_indices": [2, 5, 6]},
        ],
        "solver": {
            "type": "gmres",
            "preconditioner": {"type": "phmg"},
            "absolute_tolerance": 1.0e-7,
            "max_iterations": 10,
        },
    }

    case_file = str(tmp_path / "multigrid_profile_regions.case")
    with open(case_file, "w") as f:
        json.dump(case, f, indent=4)

    result = run_neko(launcher_script, configure_nprocs(1), case_file,
                      get_neko(), log_file)

    with open(log_file) as f:
        log = f.read()
    assert "Profile region renamed" not in log
    assert result.returncode == 0, \
        f"neko process failed with exit code {result.returncode}"

    # Each statistics row starts with the name of its region, so a row for
    # each preconditioner also shows that both were applied
    assert "Runtime statistics" in log
    stats = log.split("Runtime statistics", 1)[1]
    regions = {line.split()[0] for line in stats.splitlines() if line.strip()}
    assert "HSMG_solve" in regions
    assert "PHMG_solve" in regions
