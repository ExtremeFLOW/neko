"""Regression test for OIFS configured with scalar advection disabled."""

import json
import subprocess

import pytest

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


VALUE_TOLERANCE = {"dp": 1.0e-10, "sp": 1.0e-5}
SOLVER_TOLERANCE = {"dp": 1.0e-12, "sp": 1.0e-6}


def test_oifs_disabled_advection(
    launcher_script, log_file, tmp_path
):
    """Fluid and scalar fields retain BDF history when OIFS is inapplicable."""
    solver_tolerance = SOLVER_TOLERANCE[conftest.RP]
    mesh = tmp_path / "box.nmsh"
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
                "polynomial_order": 2,
                "dealias": False,
                "oifs": True,
            },
            "fluid": {
                "scheme": "pnpn",
                "freeze": False,
                "advection": False,
                "mu": 1.0,
                "rho": 1.0,
                "output_control": "never",
                "initial_condition": {
                    "type": "uniform",
                    "value": [1.0, 0.0, 0.0],
                },
                "velocity_solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "absolute_tolerance": solver_tolerance,
                    "max_iterations": 100,
                },
                "pressure_solver": {
                    "type": "gmres",
                    "preconditioner": {"type": "jacobi"},
                    "absolute_tolerance": solver_tolerance,
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
            "scalar": {
                "enabled": True,
                "advection": False,
                "lambda": 1.0,
                "cp": 1.0,
                "initial_condition": {"type": "uniform", "value": 1.0},
                "solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "absolute_tolerance": solver_tolerance,
                    "max_iterations": 100,
                },
                "boundary_conditions": [
                    {
                        "type": "neumann",
                        "zone_indices": [1, 2, 3, 4, 5, 6],
                        "flux": 0.0,
                    }
                ],
            },
            "simulation_components": [
                {
                    "type": "probes",
                    "compute_control": "nsamples",
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

    case_file = tmp_path / "oifs_disabled_advection.case"
    case_file.write_text(json.dumps(case, indent=2))

    result = run_neko(
        launcher_script,
        configure_nprocs(1),
        str(case_file),
        get_neko(),
        log_file,
    )
    assert result.returncode == 0

    probe_lines = (tmp_path / "probe.csv").read_text().splitlines()
    _, velocity_value, scalar_value = map(
        float, probe_lines[-1].split(",")
    )
    tolerance = VALUE_TOLERANCE[conftest.RP]
    assert velocity_value == pytest.approx(1.0, abs=tolerance)
    assert scalar_value == pytest.approx(1.0, abs=tolerance)
