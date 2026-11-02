"""Test of the GMRES space_size key of the linear solver blocks."""
from os.path import join
import json
import re

import json5
from testlib import get_neko, run_neko, configure_nprocs, get_neko_dir


def make_case(tmp_path, name, pressure, steps=2):
    """A two-step cylinder case with the given pressure solver block."""
    neko_dir = get_neko_dir()
    with open(join(neko_dir, "examples", "cylinder", "cylinder.case")) as f:
        case = json5.load(f)
    case_object = case["case"]
    case_object["mesh_file"] = join("meshes", "small_test_cyl.nmsh")
    time_object = case_object["time"]
    timestep = time_object.get("timestep", time_object.get("max_timestep"))
    time_object["end_time"] = steps * timestep
    case_object["output_directory"] = str(tmp_path)
    case_object["fluid"]["pressure_solver"].update(pressure)
    case_file = str(tmp_path / (name + ".case"))
    with open(case_file, "w") as f:
        json.dump(case, f, indent=4)
    return case_file


def run(launcher_script, case_file, log_file):
    result = run_neko(launcher_script, configure_nprocs(2), case_file,
                      get_neko(), log_file)
    assert result.returncode == 0, \
        f"neko process failed with exit code {result.returncode}"
    with open(log_file) as f:
        return f.read()


def pressure_iterations(log):
    """Pressure iterations per step from the solver summary lines."""
    return [int(m.group(1)) for m in
            re.finditer(r"\|\s*Pressure\s+(\d+)\s", log)]


def test_gmres_space_size(launcher_script, request, log_file, tmp_path):
    """space_size is read and logged, and a small space converges to
    the same tolerance as the default one, with a few more iterations."""
    base = {"type": "gmres", "preconditioner": {"type": "phmg"},
            "absolute_tolerance": 1e-6, "max_iterations": 500,
            "projection_space_size": 0}
    log_default = run(launcher_script,
                      make_case(tmp_path, request.node.name + "_30", base),
                      log_file.replace(".log", "_30.log"))
    log_small = run(launcher_script,
                    make_case(tmp_path, request.node.name + "_4",
                              dict(base, space_size=4)),
                    log_file.replace(".log", "_4.log"))
    assert "GMRES space: 30" in log_default
    assert "GMRES space: 4" in log_small
    its_default = pressure_iterations(log_default)
    its_small = pressure_iterations(log_small)
    assert its_default and its_small
    assert "did not converge" not in log_small
    # restarting costs iterations, but not many with a multigrid preconditioner
    assert its_default[-1] <= its_small[-1] <= 3 * its_default[-1]
