"""Temporal convergence of extrapolated user source terms.

A periodic box starts at rest with the scalar at zero, and a spatially uniform
user source f(t) = cos(t) is added to the x-momentum and scalar equations, with
rho = cp = 1. The fields stay uniform, so viscous, diffusive and pressure terms
vanish and both solutions reduce to u(t) = s(t) = sin(t).

The source is evaluated at the previous time levels t^n, t^{n-1}, ... and
extrapolated with the EXT scheme, so the observed convergence rate at t = 1
must match the order of the scheme. If the source is evaluated at t^{n+1}
instead, the extrapolation is shifted by one step and the error converges at
first order regardless of the scheme order.
"""

import json
import subprocess
from os.path import abspath, dirname, join

import numpy as np
import pytest

import conftest
from testlib import (
    configure_nprocs,
    get_genmeshbox,
    get_makeneko,
    run_neko,
)

USER_FILE = join(dirname(abspath(__file__)), "test_source_term_time.f90")
END_TIME = 1.0
TIMESTEPS = [0.2, 0.1, 0.05]
SOLVER_TOLERANCE = {"dp": 1.0e-12, "sp": 1.0e-6}


@pytest.fixture(scope="module")
def build(tmp_path_factory):
    """Build the user file and generate the mesh once for all orders."""
    build_dir = tmp_path_factory.mktemp("source_term_time")

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
            "2", "2", "2",
            ".true.", ".true.", ".true.",
        ],
        cwd=build_dir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    return build_dir


def make_case(mesh, output_dir, order, dt):
    solver = {
        "type": "cg",
        "preconditioner": {"type": "jacobi"},
        "absolute_tolerance": SOLVER_TOLERANCE[conftest.RP],
        "max_iterations": 200,
    }
    source = {"type": "user"}
    return {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_dir),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": END_TIME, "timestep": dt},
            "numerics": {
                "time_order": order,
                "polynomial_order": 2,
                "dealias": False,
            },
            "fluid": {
                "scheme": "pnpn",
                "advection": False,
                "mu": 1.0,
                "rho": 1.0,
                "output_control": "never",
                "initial_condition": {
                    "type": "uniform",
                    "value": [0.0, 0.0, 0.0],
                },
                "source_terms": [source],
                "velocity_solver": solver,
                "pressure_solver": solver,
            },
            "scalar": {
                "enabled": True,
                "advection": False,
                "lambda": 1.0,
                "cp": 1.0,
                "initial_condition": {"type": "uniform", "value": 0.0},
                "source_terms": [source],
                "solver": solver,
            },
            "simulation_components": [
                {
                    "type": "probes",
                    "compute_control": "tsteps",
                    "compute_value": 1,
                    "append_output": False,
                    "output_file": "probe.csv",
                    "fields": ["u", "s"],
                    "points": [
                        {
                            "type": "points",
                            "coordinates": [0.5, 0.5, 0.5],
                        }
                    ],
                }
            ],
        },
    }


@pytest.mark.parametrize("order", [2, 3])
def test_source_term_time(launcher_script, log_file, build, tmp_path, order):
    """Extrapolated user source terms converge at the order of the scheme."""
    errors = {"u": [], "s": []}

    for i, dt in enumerate(TIMESTEPS):
        run_dir = tmp_path / f"dt{i}"
        run_dir.mkdir()
        case_file = run_dir / "source_term_time.case"
        case_file.write_text(
            json.dumps(make_case(build / "box.nmsh", run_dir, order, dt),
                       indent=2)
        )

        run_log = log_file.replace(".log", f"_dt{i}.log")
        result = run_neko(
            launcher_script,
            configure_nprocs(1),
            str(case_file),
            str(build / "neko"),
            run_log,
        )
        assert result.returncode == 0, f"neko failed, see {run_log}"

        last = (run_dir / "probe.csv").read_text().splitlines()[-1]
        t, u, s = map(float, last.split(","))
        assert t == pytest.approx(END_TIME), f"last probe at t = {t}"

        exact = np.sin(t)
        errors["u"].append(abs(u - exact))
        errors["s"].append(abs(s - exact))

    for name, err in errors.items():
        rate = np.polyfit(np.log(TIMESTEPS), np.log(err), 1)[0]
        assert rate > order - 0.3, (
            f"{name}: convergence rate {rate:.2f} for time order {order}, "
            f"errors {err} for dt {TIMESTEPS}"
        )
