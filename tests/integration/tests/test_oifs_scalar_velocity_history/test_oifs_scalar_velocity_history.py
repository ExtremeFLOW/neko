"""Regression test for scalar OIFS velocity time levels."""

import json
import subprocess
from pathlib import Path

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


RESIDUAL_TOLERANCE = {"dp": 1.0e-12, "sp": 1.0e-5}
SOLVER_TOLERANCE = {"dp": 1.0e-12, "sp": 1.0e-6}


def _scalar_start_residual(log_file):
    """Return the scalar solver's initial residual from the first step."""
    for line in Path(log_file).read_text().splitlines():
        columns = line.split()
        if len(columns) >= 5 and columns[1:3] == ["|", "s"]:
            return float(columns[4])
    raise AssertionError("Scalar solver result was not found in the log")


def test_oifs_scalar_uses_matching_velocity_history(
    launcher_script, log_file, tmp_path
):
    """Scalar OIFS pairs the fluid lag field with its lagged time."""
    solver_tolerance = SOLVER_TOLERANCE[conftest.RP]
    mesh = tmp_path / "box.nmsh"
    result = subprocess.run(
        [
            get_genmeshbox(),
            "0", "1", "0", "1", "0", "1",
            "3", "3", "3",
            ".true.", ".true.", ".true.",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    assert mesh.is_file()

    case = {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(tmp_path),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": 0.1, "timestep": 0.1},
            "numerics": {
                "time_order": 1,
                "polynomial_order": 4,
                "dealias": True,
                "oifs": True,
            },
            "fluid": {
                "scheme": "pnpn",
                "rho": 1.0,
                "mu": 1.0,
                "output_control": "never",
                "initial_condition": {
                    "type": "uniform",
                    "value": [0.0, 0.0, 0.0],
                },
                "source_terms": [
                    {
                        "type": "constant",
                        "values": [1.0, 0.0, 0.0],
                    }
                ],
                "velocity_solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "projection_space_size": 0,
                    "absolute_tolerance": solver_tolerance,
                    "max_iterations": 100,
                },
                "pressure_solver": {
                    "type": "gmres",
                    "preconditioner": {"type": "hsmg"},
                    "projection_space_size": 0,
                    "absolute_tolerance": solver_tolerance,
                    "max_iterations": 800,
                },
                "boundary_conditions": [],
            },
            "scalar": {
                "enabled": True,
                "advection": True,
                "lambda": 0.0,
                "cp": 1.0,
                "initial_condition": {
                    "type": "expression",
                    "value": "sin(6.283185307179586*x)",
                },
                "solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "absolute_tolerance": solver_tolerance,
                    "max_iterations": 100,
                },
                "boundary_conditions": [],
            },
            "simulation_components": [
                {
                    "type": "probes",
                    "compute_control": "nsamples",
                    "compute_value": 1,
                    "append_output": False,
                    "output_file": "probe.csv",
                    "fields": ["u"],
                    "points": [
                        {
                            "type": "points",
                            "coordinates": [0.125, 0.25, 0.25],
                        }
                    ],
                }
            ],
        },
    }

    case_file = tmp_path / "oifs_scalar_velocity_history.case"
    case_file.write_text(json.dumps(case, indent=2))

    result = run_neko(
        launcher_script,
        configure_nprocs(1),
        str(case_file),
        get_neko(),
        log_file,
    )
    assert result.returncode == 0

    # The source accelerates the fluid during the first step. Scalar OIFS must
    # nevertheless use the initially stationary velocity stored at tlag(1).
    assert _scalar_start_residual(log_file) < RESIDUAL_TOLERANCE[conftest.RP]

    probe_lines = (tmp_path / "probe.csv").read_text().splitlines()
    _, velocity_value = map(float, probe_lines[-1].split(","))
    assert velocity_value > 5.0e-2
