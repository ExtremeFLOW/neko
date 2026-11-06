"""Regression test for runtime statistics with several scalars."""

import json
import subprocess

from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


SCALAR_NAMES = ["temperature", "dye"]


def scalar_params(name, value):
    """Parameters of a scalar with a uniform initial condition."""
    return {
        "name": name,
        "lambda": 1.0e-3,
        "cp": 1.0,
        "initial_condition": {"type": "uniform", "value": value},
        "solver": {
            "type": "cg",
            "preconditioner": {"type": "jacobi"},
            "absolute_tolerance": 1.0e-7,
            "max_iterations": 100,
        },
    }


def test_scalar_profile_regions(launcher_script, log_file, tmp_path):
    """Every scalar gets its own profiling regions in the statistics."""
    mesh = tmp_path / "box.nmsh"
    result = subprocess.run(
        [
            get_genmeshbox(),
            "0", "1", "0", "1", "0", "1",
            "2", "2", "2",
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
            "time": {"end_time": 3.0e-3, "timestep": 1.0e-3},
            "numerics": {
                "time_order": 2,
                "polynomial_order": 3,
                "dealias": False,
            },
            "runtime_statistics": {"enabled": True},
            "fluid": {
                "scheme": "pnpn",
                "freeze": True,
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
                    "absolute_tolerance": 1.0e-7,
                    "max_iterations": 100,
                },
                "pressure_solver": {
                    "type": "gmres",
                    "preconditioner": {"type": "jacobi"},
                    "absolute_tolerance": 1.0e-7,
                    "max_iterations": 100,
                },
            },
            "scalars": [
                scalar_params(name, float(i + 1))
                for i, name in enumerate(SCALAR_NAMES)
            ],
        },
    }
    case_file = tmp_path / "scalar_profile_regions.case"
    case_file.write_text(json.dumps(case, indent=2))

    result = run_neko(
        launcher_script,
        configure_nprocs(1),
        str(case_file),
        get_neko(),
        log_file,
    )
    assert result.returncode == 0

    # Each statistics row starts with the name of its region
    with open(log_file) as f:
        log = f.read()
    assert "Runtime statistics" in log
    stats = log.split("Runtime statistics", 1)[1]
    regions = {line.split()[0] for line in stats.splitlines() if line.strip()}
    for name in SCALAR_NAMES:
        for suffix in ["", "_residual", "_solve"]:
            assert name + suffix in regions
