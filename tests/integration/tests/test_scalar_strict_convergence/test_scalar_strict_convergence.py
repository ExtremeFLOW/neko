"""Tests of the strict_convergence key of a scalar.

Each scalar solve is capped at one iteration, far from its tolerance, so it
never converges. Without the key the run only logs that. With the key set on
a scalar, the run must stop at the first solve of that scalar, and only that
scalar's key counts.
"""

import json
import subprocess

import pytest

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


SOLVER_TOLERANCE = {"dp": 1.0e-10, "sp": 1.0e-5}


@pytest.fixture
def mesh(tmp_path):
    """A box of 2x2x2 elements, with its six faces in zones 1 to 6."""
    result = subprocess.run(
        [
            get_genmeshbox(),
            "0", "1", "0", "1", "0", "1",
            "2", "2", "2",
            ".false.", ".false.", ".false.",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    mesh_file = tmp_path / "box.nmsh"
    assert mesh_file.is_file()
    return mesh_file


def scalar(name, max_iterations, strict_convergence=None):
    """A diffusing scalar, driven by its value on the inlet face."""
    scalar_object = {
        "name": name,
        "lambda": 1.0,
        "cp": 1.0,
        "initial_condition": {"type": "uniform", "value": 0.0},
        "solver": {
            "type": "cg",
            "preconditioner": {"type": "jacobi"},
            "absolute_tolerance": SOLVER_TOLERANCE[conftest.RP],
            "max_iterations": max_iterations,
        },
        "boundary_conditions": [
            {"type": "dirichlet", "value": 1.0, "zone_indices": [1]},
            {
                "type": "neumann",
                "flux": 0.0,
                "zone_indices": [2, 3, 4, 5, 6],
            },
        ],
    }
    if strict_convergence is not None:
        scalar_object["strict_convergence"] = strict_convergence
    return scalar_object


def run_case(tmp_path, mesh, name, scalars, launcher_script, log_file):
    """Run three steps with the fluid frozen, and return the result and log.

    A single scalar goes in `case.scalar`, several in `case.scalars`.
    """
    case = {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(tmp_path),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": 0.03, "timestep": 0.01},
            "numerics": {
                "time_order": 1,
                "polynomial_order": 3,
                "dealias": False,
            },
            "fluid": {
                "scheme": "pnpn",
                "freeze": True,
                "Re": 1.0,
                "output_control": "never",
                "initial_condition": {
                    "type": "uniform",
                    "value": [1.0, 0.0, 0.0],
                },
                "velocity_solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "absolute_tolerance": 1.0e-6,
                    "max_iterations": 100,
                },
                "pressure_solver": {
                    "type": "gmres",
                    "preconditioner": {"type": "jacobi"},
                    "absolute_tolerance": 1.0e-6,
                    "max_iterations": 100,
                },
                "boundary_conditions": [
                    {
                        "type": "velocity_value",
                        "zone_indices": [1, 2, 3, 4, 5, 6],
                        "value": [1.0, 0.0, 0.0],
                    }
                ],
            },
        },
    }
    if len(scalars) == 1:
        case["case"]["scalar"] = scalars[0]
    else:
        case["case"]["scalars"] = scalars

    case_file = tmp_path / f"{name}.case"
    case_file.write_text(json.dumps(case, indent=2))

    run_log = log_file.replace(".log", f"_{name}.log")
    result = run_neko(
        launcher_script, configure_nprocs(1), str(case_file), get_neko(),
        run_log
    )
    with open(run_log) as f:
        return result, f.read()


def test_scalar_not_strict(launcher_script, log_file, tmp_path, mesh):
    """By default, a scalar solve that does not converge is only logged."""
    result, log = run_case(
        tmp_path, mesh, "default", [scalar("tracer", 1)], launcher_script,
        log_file
    )
    assert result.returncode == 0, \
        f"neko process failed with exit code {result.returncode}"
    assert "Scalar solver did not converge for tracer" in log
    assert "ERROR" not in log


def test_scalar_strict(launcher_script, log_file, tmp_path, mesh):
    """With strict_convergence, the same run stops with an error."""
    result, log = run_case(
        tmp_path, mesh, "strict", [scalar("tracer", 1, True)],
        launcher_script, log_file
    )
    assert result.returncode != 0, \
        "neko did not stop on a scalar solve that did not converge"
    assert "ERROR: Scalar solver did not converge for tracer" in log


def test_scalars_strict_per_scalar(launcher_script, log_file, tmp_path, mesh):
    """In `case.scalars`, each scalar's own key decides for that scalar."""
    # A strict scalar that converges does not stop the run, and the scalar
    # that does not converge is not strict.
    result, log = run_case(
        tmp_path, mesh, "other_not_strict",
        [scalar("converged", 800, True), scalar("capped", 1)],
        launcher_script, log_file
    )
    assert result.returncode == 0, \
        f"neko process failed with exit code {result.returncode}"
    assert "Scalar solver did not converge for capped" in log
    assert "did not converge for converged" not in log
    assert "ERROR" not in log

    # Setting the key on the scalar that does not converge stops the run.
    result, log = run_case(
        tmp_path, mesh, "capped_strict",
        [scalar("converged", 800), scalar("capped", 1, True)],
        launcher_script, log_file
    )
    assert result.returncode != 0, \
        "neko did not stop on a scalar solve that did not converge"
    assert "ERROR: Scalar solver did not converge for capped" in log
