"""Integration coverage for scaling the initial velocity to a forced flow rate.

With ``flow_rate_force`` the initial velocity is scaled to the target flow
rate before the first step, so that the forcing does not have to correct it
impulsively. A field without flow in the forced direction is left to the
forcing, and one already at the target is left alone. The runs take one step
in a fully periodic box, where a uniform field is a steady solution, so the
probes see the scaled initial condition, or what the forcing made of it.
"""

import json
import subprocess
from pathlib import Path

import pytest

import conftest
from testlib import (
    configure_nprocs,
    get_genmeshbox,
    get_neko,
    run_neko,
)

NPROCS = 2
ELEMENTS = (3, 3, 3)
BULK_VELOCITY = 1.0

# Tolerance of the probed values, and of the linear solves of the one step.
TOLERANCE = {"dp": 1.0e-6, "sp": 1.0e-4}
SOLVER_TOLERANCE = {"dp": 1.0e-10, "sp": 1.0e-6}
POINTS = [(0.37, 0.61, 0.25), (0.8, 0.2, 0.7)]


def _solver(solver_type):
    """Return a Jacobi preconditioned solver configuration."""
    return {
        "type": solver_type,
        "preconditioner": {"type": "jacobi"},
        "absolute_tolerance": SOLVER_TOLERANCE[conftest.RP],
        "max_iterations": 500,
    }


def _case(mesh, output_directory, initial_velocity):
    """Build a one-step, fully periodic box case with a forced flow rate."""
    return {
        "version": 1.0,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_directory),
            "output_boundary": False,
            "output_checkpoints": False,
            "output_at_end": False,
            "time": {"end_time": 1.0e-3, "timestep": 1.0e-3},
            "numerics": {
                "time_order": 3,
                "polynomial_order": 5,
                "dealias": False,
            },
            "fluid": {
                "scheme": "pnpn",
                "Re": 100.0,
                "initial_condition": {
                    "type": "uniform",
                    "value": initial_velocity,
                },
                "flow_rate_force": {
                    "direction": 1,
                    "value": BULK_VELOCITY,
                    "use_averaged_flow": True,
                },
                "velocity_solver": _solver("cg"),
                # The pressure problem is pure Neumann; cg drifts on it in
                # single precision, gmres does not.
                "pressure_solver": _solver("gmres"),
                "output_control": "never",
            },
            "simulation_components": [
                {
                    "type": "probes",
                    "compute_control": "nsamples",
                    "compute_value": 1,
                    "append_output": False,
                    "output_file": "probe.csv",
                    "fields": ["u", "v", "w"],
                    "points": [{
                        "type": "points",
                        "coordinates": [c for p in POINTS for c in p],
                    }],
                }
            ],
        },
    }


def _read_probes(probe_file, npoints):
    """Return the last sampled (u, v, w) of every probe."""
    rows = [
        line.strip()
        for line in probe_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    samples = rows[1 + npoints:]
    assert len(samples) >= npoints, rows
    return [[float(e) for e in row.split(",")][1:]
            for row in samples[-npoints:]]


def _run(assets, initial_velocity, name):
    """Run one case and return the probed values and the log."""
    run_dir = assets["workdir"] / name
    run_dir.mkdir()
    case_file = run_dir / f"{name}.case"
    case_file.write_text(
        json.dumps(_case(assets["mesh"], run_dir, initial_velocity), indent=2)
        + "\n",
        encoding="utf-8",
    )
    log_file = run_dir / f"{name}.log"
    result = run_neko(
        str(assets["launcher"]),
        configure_nprocs(NPROCS),
        str(case_file),
        str(assets["neko"]),
        str(log_file),
    )
    log = log_file.read_text(encoding="utf-8", errors="replace")
    assert result.returncode == 0, log
    return _read_probes(run_dir / "probe.csv", len(POINTS)), log


@pytest.fixture(scope="module")
def flow_rate_assets(tmp_path_factory, request):
    """Generate the fully periodic unit cube shared by the runs."""
    workdir = tmp_path_factory.mktemp("flow_rate_ic")
    genmeshbox = Path(get_genmeshbox()).resolve()
    result = subprocess.run(
        [
            str(genmeshbox),
            "0", "1", "0", "1", "0", "1",
            *(str(value) for value in ELEMENTS),
            ".true.", ".true.", ".true.",
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    return {
        "workdir": workdir,
        "neko": Path(get_neko()).resolve(),
        "launcher": Path(
            request.config.getoption("--launcher-script")).resolve(),
        "mesh": workdir / "box.nmsh",
    }


def test_initial_condition_is_scaled_to_the_flow_rate(flow_rate_assets):
    """A uniform field at half the bulk velocity is scaled up to it."""
    tolerance = TOLERANCE[conftest.RP]
    values, log = _run(flow_rate_assets, [0.5, 0.2, 0.0], "scaled")
    assert "scaled to" in log
    for u, v, w in values:
        assert u == pytest.approx(BULK_VELOCITY, abs=tolerance)
        assert v == pytest.approx(0.4, abs=tolerance)
        assert w == pytest.approx(0.0, abs=tolerance)


def test_initial_condition_at_rest_is_left_to_the_forcing(flow_rate_assets):
    """A field without flow in the forced direction is not scaled; the forcing
    brings it to the bulk velocity in the first step instead."""
    tolerance = TOLERANCE[conftest.RP]
    values, log = _run(flow_rate_assets, [0.0, 0.3, 0.0], "at_rest")
    assert "not scaled" in log
    assert "experimental" not in log
    for u, v, w in values:
        assert u == pytest.approx(BULK_VELOCITY, abs=tolerance)
        assert v == pytest.approx(0.3, abs=tolerance)
        assert w == pytest.approx(0.0, abs=tolerance)


def test_initial_condition_at_the_flow_rate_is_left_alone(flow_rate_assets):
    """A field at the bulk velocity is left alone, without a warning."""
    tolerance = TOLERANCE[conftest.RP]
    values, log = _run(flow_rate_assets, [BULK_VELOCITY, 0.3, 0.0],
                       "at_target")
    assert "at the target" in log
    assert "experimental" not in log
    for u, v, w in values:
        assert u == pytest.approx(BULK_VELOCITY, abs=tolerance)
        assert v == pytest.approx(0.3, abs=tolerance)
        assert w == pytest.approx(0.0, abs=tolerance)
