from os.path import join
from testlib import get_genmeshbox, get_makeneko, run_neko, configure_nprocs
from conftest import RP
import subprocess
import json
import re

# Absolute tolerances of the linear solvers and of the check on s(dt).
SOLVER_TOLERANCE = {"dp": 1e-12, "sp": 1e-6}
CHECK_TOLERANCE = {"dp": 1e-10, "sp": 1e-5}


def test_scalar_implicit_source(launcher_script, log_file, tmp_path):
    """One BDF1 step of rho cp ds/dt = rho cp (a s + b) with a = -90, b = 2,
    s(0) = 1, dt = 0.1, on a periodic box with a uniform scalar. With a s
    treated implicitly, s(dt) = (1 + b dt) / (1 - a dt) = 0.12, independent
    of rho cp. An explicit treatment would give 1 + (a + b) dt = -7.8."""

    makeneko = get_makeneko()
    result = subprocess.run(
        [makeneko, join("tests", "test_scalar_implicit_source",
                        "test_scalar_implicit_source.f90")],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True)
    assert result.returncode == 0, result.stdout

    genmeshbox = get_genmeshbox()
    result = subprocess.run(
        [genmeshbox, "0", "1", "0", "1", "0", "1", "3", "3", "3",
         ".true.", ".true.", ".true."],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True)
    assert result.returncode == 0, result.stdout

    case_template = join("tests", "test_scalar_implicit_source",
                         "test_scalar_implicit_source.case")
    with open(case_template, "r", encoding="utf-8") as handle:
        case_data = json.load(handle)
    case_data["case"]["output_directory"] = str(tmp_path)
    for solver in (case_data["case"]["fluid"]["velocity_solver"],
                   case_data["case"]["fluid"]["pressure_solver"],
                   case_data["case"]["scalar"]["solver"]):
        solver["absolute_tolerance"] = SOLVER_TOLERANCE[RP]

    case_file = tmp_path / "test_scalar_implicit_source.case"
    with open(case_file, "w", encoding="utf-8") as handle:
        json.dump(case_data, handle, indent=2)
        handle.write("\n")

    result = run_neko(launcher_script, configure_nprocs(1), case_file,
                      "./neko", log_file)
    assert result.returncode == 0, \
        f"neko process failed with exit code {result.returncode}"

    with open(log_file, "r", encoding="utf-8") as handle:
        log = handle.read()
    values = re.findall(
        r"implicit_source_scalar\s+(\S+)\s+(\S+)", log)
    assert values, "Scalar value was not reported"

    s_min, s_max = (float(v) for v in values[-1])
    assert abs(s_min - 0.12) < CHECK_TOLERANCE[RP], f"min(s) = {s_min}"
    assert abs(s_max - 0.12) < CHECK_TOLERANCE[RP], f"max(s) = {s_max}"
