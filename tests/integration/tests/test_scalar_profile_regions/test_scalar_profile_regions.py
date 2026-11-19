"""Regression test for runtime statistics with several scalars."""

import json
import re
import subprocess

from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


SCALAR_NAMES = ["temperature", "dye"]

# Regions shared by all scalars, with fixed ids like the fluid's
COMMON_REGIONS = ["Scalar", "Scalar_residual", "Scalar_solve"]

# A statistics row: the region name and three times
STATS_ROW = re.compile(
    r"^\s*(\S.*?)\s+(\S+E[+-]\d+)\s+(\S+E[+-]\d+)\s+(\S+E[+-]\d+)\s*$"
)


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
    """The scalars share common profiling regions and each has its own."""
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

    with open(log_file) as f:
        log = f.read()
    assert "Normal end." in log
    assert "Runtime statistics" in log
    stats = log.split("Runtime statistics", 1)[1]
    regions = [
        match.group(1).strip()
        for match in map(STATS_ROW.match, stats.splitlines())
        if match
    ]

    # One row for each common region, whatever the number of scalars
    for name in COMMON_REGIONS:
        assert regions.count(name) == 1, regions

    # One row for each region of every scalar
    own_regions = [
        name + suffix
        for name in SCALAR_NAMES
        for suffix in ["", "_residual", "_solve"]
    ]
    for name in own_regions:
        assert regions.count(name) == 1, regions

    # The report lists the regions by id: the fixed ids first, then those
    # allocated by name in the order they are first started. "Time-Step"
    # has no id and starts before any scalar region, so a common region
    # listed after it has no fixed id.
    time_step = regions.index("Time-Step")
    for name in COMMON_REGIONS:
        assert regions.index(name) < time_step, regions
