"""Regression test for checkpoint restarts with OIFS advection."""

import json
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


SOLVER_TOLERANCE = {"dp": 1.0e-12, "sp": 1.0e-6}
TWO_PI = "6.283185307179586"


def _case(mesh, output_directory, target_cfl):
    """Return a periodic Taylor-Green vortex case with a passive scalar."""
    tolerance = SOLVER_TOLERANCE[conftest.RP]
    return {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_directory),
            "output_at_end": True,
            "output_boundary": False,
            "output_checkpoints": True,
            "checkpoint_control": "tsteps",
            "checkpoint_value": 10,
            # dt gives a CFL number of about 1.15 on this mesh.
            "time": {"end_time": 0.8, "timestep": 0.04},
            "numerics": {
                "time_order": 3,
                "polynomial_order": 3,
                "dealias": True,
                "oifs": True,
                "oifs_target_cfl": target_cfl,
            },
            "fluid": {
                "scheme": "pnpn",
                "rho": 1.0,
                "mu": 1.0e-3,
                "output_control": "never",
                "initial_condition": {
                    "type": "expression",
                    "value": [
                        f"sin({TWO_PI}*x)*cos({TWO_PI}*y)",
                        f"-cos({TWO_PI}*x)*sin({TWO_PI}*y)",
                        "0",
                    ],
                },
                "velocity_solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "projection_space_size": 0,
                    "absolute_tolerance": tolerance,
                    "max_iterations": 200,
                },
                # Jacobi keeps the device runs reproducible, see the note on
                # TreeAMG in test_scalar_restart.
                "pressure_solver": {
                    "type": "gmres",
                    "preconditioner": {"type": "jacobi"},
                    "projection_space_size": 0,
                    "absolute_tolerance": tolerance,
                    "max_iterations": 800,
                },
                "boundary_conditions": [],
            },
            "scalar": {
                "enabled": True,
                "advection": True,
                "lambda": 1.0e-3,
                "cp": 1.0,
                "initial_condition": {
                    "type": "expression",
                    "value": f"sin({TWO_PI}*x)",
                },
                "solver": {
                    "type": "cg",
                    "preconditioner": {"type": "jacobi"},
                    "absolute_tolerance": tolerance,
                    "max_iterations": 200,
                },
                "boundary_conditions": [],
            },
        },
    }


def _run(launcher_script, case, case_file, log_file):
    """Write the case file and run Neko on it."""
    case_file.write_text(json.dumps(case, indent=2))
    result = run_neko(
        launcher_script,
        configure_nprocs(1),
        str(case_file),
        get_neko(),
        str(log_file),
    )
    assert result.returncode == 0, f"Neko failed, see {log_file}"


# Target CFL 2 gives 1 RK4 substep per time step, target CFL 4 gives 2.
@pytest.mark.parametrize("target_cfl", [2.0, 4.0])
def test_oifs_restart(
    launcher_script, request, tmp_path, monkeypatch, target_cfl
):
    """A restarted OIFS run reproduces the uninterrupted run bit for bit."""
    # On devices the autotuner picks the kernel formulation per run, and the
    # 1D and KSTEP Ax kernels round differently. Pin the formulation so that
    # both runs use the same kernels. The CPU backend ignores this variable.
    monkeypatch.setenv("NEKO_AUTOTUNE", "1D")

    result = subprocess.run(
        [
            get_genmeshbox(),
            "0", "1", "0", "1", "0", "0.125",
            "8", "8", "1",
            ".true.", ".true.", ".true.",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    mesh = tmp_path / "box.nmsh"

    # Keep the restart path short enough for Neko's fixed-size log buffer.
    full_dir_context = TemporaryDirectory(prefix="neko-oifs-", dir="/tmp")
    request.addfinalizer(full_dir_context.cleanup)
    full_dir = Path(full_dir_context.name)
    restart_dir_context = TemporaryDirectory(prefix="neko-oifs-", dir="/tmp")
    request.addfinalizer(restart_dir_context.cleanup)
    restart_dir = Path(restart_dir_context.name)

    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    name = request.node.name

    # Uninterrupted run, writing a checkpoint after 10 and 20 steps.
    case = _case(mesh, full_dir, target_cfl)
    _run(
        launcher_script,
        case,
        tmp_path / "full.case",
        log_dir / f"{name}_full.log",
    )

    # Restart from the checkpoint after 10 steps and run to the same end time.
    case = _case(mesh, restart_dir, target_cfl)
    case["case"]["restart_file"] = str(full_dir / "fluid00000.chkp")
    case["case"]["output_checkpoints"] = False
    del case["case"]["checkpoint_control"]
    del case["case"]["checkpoint_value"]
    _run(
        launcher_script,
        case,
        tmp_path / "restart.case",
        log_dir / f"{name}_restart.log",
    )

    # The field files hold u, v, w, p and the scalar at the end time.
    full = (full_dir / "field0.f00000").read_bytes()
    restarted = (restart_dir / "field0.f00000").read_bytes()
    assert restarted == full, (
        "The restarted run does not reproduce the uninterrupted run"
    )
