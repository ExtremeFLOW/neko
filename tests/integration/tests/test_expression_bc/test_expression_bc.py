"""Integration coverage for the expression boundary conditions.

``expression_velocity`` and the scalar ``expression_dirichlet`` are strong
Dirichlet conditions, so the value they prescribe has to survive the linear
solve, not merely be written into the field before it. A condition that is
applied but not registered as Dirichlet with the scheme is left unconstrained
by the solver and drifts away within one step, silently, which is what this
test guards against.
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

# The prescribed boundary values.
INFLOW_U = 1.0
SCALAR_VALUE = 3.5

TOLERANCE = {"dp": 1.0e-9, "sp": 1.0e-5}


def _solver(solver_type):
    return {
        "type": solver_type,
        "preconditioner": {"type": "jacobi"},
        "absolute_tolerance": 1.0e-10,
        "max_iterations": 500,
    }


def _case(mesh, output_directory, expression):
    """Build a duct case, prescribing the inflow either by value or by
    expression. Both spell out exactly the same boundary values."""
    if expression:
        inflow = {
            "type": "expression_velocity",
            "zone_indices": [1],
            "value": [str(INFLOW_U), "0", "0"],
        }
        scalar_inflow = {
            "type": "expression_dirichlet",
            "zone_indices": [1],
            "value": str(SCALAR_VALUE),
        }
    else:
        inflow = {
            "type": "velocity_value",
            "zone_indices": [1],
            "value": [INFLOW_U, 0.0, 0.0],
        }
        scalar_inflow = {
            "type": "dirichlet",
            "zone_indices": [1],
            "value": SCALAR_VALUE,
        }

    return {
        "version": 1.0,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_directory),
            "output_boundary": False,
            "output_checkpoints": False,
            "output_at_end": False,
            "time": {"end_time": 1.0e-2, "timestep": 1.0e-3},
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
                    "value": [INFLOW_U, 0.0, 0.0],
                },
                "boundary_conditions": [
                    inflow,
                    {"type": "outflow", "zone_indices": [2]},
                    {"type": "no_slip", "zone_indices": [3, 4]},
                ],
                "velocity_solver": _solver("cg"),
                "pressure_solver": _solver("gmres"),
                "output_control": "never",
            },
            "scalar": {
                "enabled": True,
                "Pe": 1.0,
                "initial_condition": {"type": "uniform", "value": 0.0},
                "boundary_conditions": [
                    scalar_inflow,
                    {"type": "neumann", "zone_indices": [2, 3, 4],
                     "flux": 0.0},
                ],
                "solver": _solver("cg"),
            },
            "simulation_components": [
                {
                    "type": "probes",
                    "compute_control": "tsteps",
                    "compute_value": 1,
                    "append_output": False,
                    "output_file": "probe.csv",
                    "fields": ["u", "s"],
                    # On the inflow plane, where both conditions are imposed.
                    "points": [{"type": "points",
                                "coordinates": [0.0, 0.5, 0.5]}],
                }
            ],
        },
    }


def _generate_mesh(genmeshbox, workdir):
    """Generate a unit cube, periodic in z only."""
    result = subprocess.run(
        [
            str(genmeshbox),
            "0", "1", "0", "1", "0", "1",
            *(str(value) for value in ELEMENTS),
            ".false.", ".false.", ".true.",
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    return workdir / "box.nmsh"


def _run(assets, expression, name):
    """Run one case and return every sampled (u, s) pair."""
    run_dir = assets["workdir"] / name
    run_dir.mkdir()
    case_file = run_dir / f"{name}.case"
    case_file.write_text(
        json.dumps(_case(assets["mesh"], run_dir, expression), indent=2)
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

    rows = [
        line.strip()
        for line in (run_dir / "probe.csv").read_text(
            encoding="utf-8").splitlines()
        if line.strip()
    ]
    # One header line, one line of coordinates, then time, u, s per sample.
    samples = [[float(e) for e in r.split(",")] for r in rows[2:]]
    assert len(samples) >= 5, rows
    return [(s[1], s[2]) for s in samples]


@pytest.fixture(scope="module")
def expression_bc_assets(tmp_path_factory, request):
    workdir = tmp_path_factory.mktemp("expression_bc")
    neko = Path(get_neko()).resolve()
    genmeshbox = Path(get_genmeshbox()).resolve()
    launcher = Path(request.config.getoption("--launcher-script")).resolve()
    return {
        "workdir": workdir,
        "neko": neko,
        "launcher": launcher,
        "mesh": _generate_mesh(genmeshbox, workdir),
    }


def test_expression_boundary_conditions_are_enforced(expression_bc_assets):
    """An expression Dirichlet value has to hold for every step, and match the
    plain Dirichlet condition prescribing the same thing."""
    tolerance = TOLERANCE[conftest.RP]
    by_value = _run(expression_bc_assets, False, "by_value")
    by_expression = _run(expression_bc_assets, True, "by_expression")

    for step, ((u_v, s_v), (u_e, s_e)) in enumerate(
            zip(by_value, by_expression), start=1):
        assert u_v == pytest.approx(INFLOW_U, abs=tolerance), f"step {step}"
        assert s_v == pytest.approx(SCALAR_VALUE, abs=tolerance), f"step {step}"
        assert u_e == pytest.approx(INFLOW_U, abs=tolerance), f"step {step}"
        assert s_e == pytest.approx(SCALAR_VALUE, abs=tolerance), f"step {step}"
