"""Units of the weak tractions in the fluid.

A box is periodic in x and z. On the boundaries y = 0 and y = 1, the velocity
has a Neumann condition with the traction [tau, 0, 0] and the pressure is zero.
The fluid starts at rest with no other forcing and without advection, so
u = u(y) and v = w = p = 0. The viscous term integrates to the boundary
traction, so the x-momentum balance gives

    rho * d<u>/dt * V = 2 * tau * A,

with A = V = 1, i.e. <u>(t) = 2 tau t / rho. The traction is a stress, so the
volume average must scale with 1 / rho.
"""

import json
import re
import subprocess
from os.path import abspath, dirname, join

import pytest

import conftest
from testlib import (
    configure_nprocs,
    get_genmeshbox,
    get_makeneko,
    run_neko,
)

USER_FILE = join(dirname(abspath(__file__)), "test_fluid_traction_rho.f90")
END_TIME = 0.1
TAU = 1.0
SOLVER_TOLERANCE = {"dp": 1.0e-12, "sp": 1.0e-6}
RELATIVE_TOLERANCE = {"dp": 1.0e-8, "sp": 1.0e-4}


@pytest.fixture(scope="module")
def build(tmp_path_factory):
    """Build the user file and generate the mesh once for all densities."""
    build_dir = tmp_path_factory.mktemp("fluid_traction_rho")

    result = subprocess.run(
        [get_makeneko(), USER_FILE],
        cwd=build_dir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    result = subprocess.run(
        [
            get_genmeshbox(),
            "0", "1", "0", "1", "0", "1",
            "3", "3", "3",
            ".true.", ".false.", ".true.",
        ],
        cwd=build_dir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    return build_dir


def make_case(mesh, output_dir, rho):
    def solver(type, preconditioner):
        return {
            "type": type,
            "preconditioner": {"type": preconditioner},
            "projection_space_size": 0,
            "absolute_tolerance": SOLVER_TOLERANCE[conftest.RP],
            "max_iterations": 500,
        }

    return {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_dir),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": END_TIME, "timestep": 0.01},
            "numerics": {
                "time_order": 2,
                "polynomial_order": 3,
                "dealias": False,
            },
            "fluid": {
                "scheme": "pnpn",
                "advection": False,
                "mu": 1.0,
                "rho": rho,
                "output_control": "never",
                "initial_condition": {
                    "type": "uniform",
                    "value": [0.0, 0.0, 0.0],
                },
                "velocity_solver": solver("cg", "jacobi"),
                "pressure_solver": solver("gmres", "hsmg"),
                "boundary_conditions": [
                    {
                        "velocity": {
                            "type": "neumann",
                            "flux": [TAU, 0.0, 0.0],
                        },
                        "pressure": {"type": "zero_dirichlet"},
                        "zone_indices": [3, 4],
                    }
                ],
            },
        },
    }


def mean_velocity(launcher_script, log_file, build, run_dir, rho):
    """Run the case with density rho and return the final volume-averaged u."""
    run_dir.mkdir()
    case_file = run_dir / "fluid_traction_rho.case"
    case_file.write_text(
        json.dumps(make_case(build / "box.nmsh", run_dir, rho), indent=2)
    )

    result = run_neko(
        launcher_script,
        configure_nprocs(1),
        str(case_file),
        str(build / "neko"),
        log_file,
    )
    assert result.returncode == 0, f"neko failed, see {log_file}"

    with open(log_file) as f:
        matches = re.findall(r"Mean u:\s+(\S+)\s+(\S+)", f.read())
    assert matches, f"no mean velocity in {log_file}"
    t, u_mean = map(float, matches[-1])
    assert t == pytest.approx(END_TIME), f"last mean velocity at t = {t}"
    return u_mean


@pytest.mark.parametrize("rho", [1.0, 2.0])
def test_fluid_traction_rho(launcher_script, log_file, build, tmp_path, rho):
    """The volume-averaged velocity driven by a boundary traction is
    2 tau t / rho."""
    u_mean = mean_velocity(launcher_script, log_file, build,
                           tmp_path / "run", rho)

    expected = 2.0 * TAU * END_TIME / rho
    assert u_mean == pytest.approx(
        expected, rel=RELATIVE_TOLERANCE[conftest.RP]
    ), f"mean u {u_mean} at t = {END_TIME}, expected {expected} for rho = {rho}"
