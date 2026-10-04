"""Tests of the linear-solver keys and start-up checks: the GMRES space size,
the operator and preconditioner check with its compatibility verdict, and the
true-residual check."""
from os.path import join
import json
import re

import json5
import pytest
from testlib import get_neko, run_neko, configure_nprocs, get_neko_dir


def make_case(tmp_path, name, pressure, velocity=None, steps=2):
    """A two-step cylinder case with the given solver blocks."""
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
    if velocity:
        case_object["fluid"]["velocity_solver"].update(velocity)
    case_file = join("tests", "test_linear_solvers", name + ".case")
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
    """gmres_space_size is read and logged, and a small space converges to
    the same tolerance as the default one, with a few more iterations."""
    base = {"type": "gmres", "preconditioner": {"type": "phmg"},
            "absolute_tolerance": 1e-6, "max_iterations": 500,
            "projection_space_size": 0}
    log_default = run(launcher_script,
                      make_case(tmp_path, request.node.name + "_30", base),
                      log_file.replace(".log", "_30.log"))
    log_small = run(launcher_script,
                    make_case(tmp_path, request.node.name + "_4",
                              dict(base, gmres_space_size=4)),
                    log_file.replace(".log", "_4.log"))
    assert "GMRES space: 30" in log_default
    assert "GMRES space: 4" in log_small
    its_default = pressure_iterations(log_default)
    its_small = pressure_iterations(log_small)
    assert its_default and its_small
    assert "did not converge" not in log_small
    # restarting costs iterations, but not many with a multigrid preconditioner
    assert its_default[-1] <= its_small[-1] <= 3 * its_default[-1]


def test_solver_check_compatible(launcher_script, request, log_file, tmp_path):
    """cg with the default phmg passes the start-up check, pressure and
    velocity blocks are printed, and the true residual matches."""
    log = run(launcher_script,
              make_case(tmp_path, request.node.name,
                        {"type": "cg", "preconditioner": {"type": "phmg"},
                         "absolute_tolerance": 1e-6,
                         "residual_check_interval": 1}),
              log_file)
    assert "Solver check: Pressure" in log
    assert "Solver check: Velocity" in log
    assert "cg + phmg is compatible" in log
    assert "Pressure true residual" in log
    assert "does not measure the true residual" not in log
    assert "WARNING" not in log.split("Solver check: Pressure")[1] \
        .split("Solver check: Velocity")[0]


def test_solver_check_incompatible(launcher_script, request, log_file,
                                   tmp_path):
    """cg with hsmg, whose coarse solve is a Krylov iteration, is reported as
    not compatible and the run continues."""
    log = run(launcher_script,
              make_case(tmp_path, request.node.name,
                        {"type": "cg", "preconditioner": {"type": "hsmg"},
                         "absolute_tolerance": 1e-5, "max_iterations": 2000}),
              log_file)
    assert "Solver check: Pressure" in log
    assert "requires" in log and "hsmg is not here" in log
    assert "Normal end." in log


def test_solver_check_gmres_accepts_all(launcher_script, request, log_file,
                                        tmp_path):
    """gmres has no requirements: hsmg is reported compatible with it."""
    log = run(launcher_script,
              make_case(tmp_path, request.node.name,
                        {"type": "gmres", "preconditioner": {"type": "hsmg"},
                         "absolute_tolerance": 1e-5}),
              log_file)
    assert "gmres + hsmg is compatible" in log
