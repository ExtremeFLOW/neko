"""Regression test for restarts that interpolate the checkpoint.

A checkpoint is interpolated when it is read on another mesh, given by
`case.restart_mesh_file`, or at another polynomial order. The extrapolated
right-hand-side history in the checkpoint (abx1, abx2, ...) is weighted by the
mass matrix of the discretisation it was written on, so it must not be carried
over. Before this was handled, the velocity blew up within a single time step
after a restart on a finer mesh.
"""

import json
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pytest

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_neko, run_neko


SOLVER_TOLERANCE = {"dp": 1.0e-10, "sp": 1.0e-6}
TWO_PI = "6.283185307179586"
TIMESTEP = 0.01
# The checkpoint is written after 2 steps, the restart then takes 2 more.
CHECKPOINT_STEPS = 2
END_TIME = 4 * TIMESTEP
# Largest allowed difference between the restarted run and a run on the
# target discretisation from t = 0. Both fields have unit amplitude. With the
# fix the difference is about 2e-3; before it, it was about 40 for a change of
# mesh and 6e-2 for a change of polynomial order.
MAX_DIFFERENCE = 1.0e-2


def _case(mesh, output_directory, order):
    """Return a periodic Taylor-Green vortex case with a passive scalar."""
    tolerance = SOLVER_TOLERANCE[conftest.RP]
    return {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_directory),
            "output_at_end": True,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": END_TIME, "timestep": TIMESTEP},
            "numerics": {
                "time_order": 3,
                "polynomial_order": order,
                "dealias": True,
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


def _genmeshbox(directory, nel, name):
    """Write a periodic unit box with nel^3 elements to directory/name."""
    result = subprocess.run(
        [
            get_genmeshbox(),
            "0", "1", "0", "1", "0", "1",
            str(nel), str(nel), str(nel),
            ".true.", ".true.", ".true.",
        ],
        cwd=directory,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    mesh = directory / name
    (directory / "box.nmsh").rename(mesh)
    return mesh


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


def _read_fld(path):
    """Return the velocity and scalar of a 3D fld file written with the mesh.

    The data follows the 132 byte header, a 4 byte endianness marker and the
    element numbering. Each record is stored element by element, the vector
    ones with all x values of an element before its y and z values.
    """
    data = path.read_bytes()
    header = data[:132].decode().split()
    wdsz = int(header[1])
    lxyz = int(header[2]) * int(header[3]) * int(header[4])
    nel = int(header[5])
    rdcode = header[11]
    assert rdcode.startswith("XUP"), f"Unexpected fields {rdcode} in {path}"

    dtype = np.float32 if wdsz == 4 else np.float64
    offset = 132 + 4 + 4 * nel
    values = np.frombuffer(data, dtype=dtype, offset=offset)

    fields = {}
    n_vector = 3 * lxyz * nel
    # Skip the coordinates
    start = n_vector
    velocity = values[start:start + n_vector].reshape(nel, 3, lxyz)
    fields["u"] = velocity[:, 0, :]
    fields["v"] = velocity[:, 1, :]
    fields["w"] = velocity[:, 2, :]
    # Skip the pressure, which is only determined up to a constant
    start += n_vector + lxyz * nel
    if "T" in rdcode:
        fields["s"] = values[start:start + lxyz * nel].reshape(nel, lxyz)
    return fields


@pytest.mark.parametrize("change", ["mesh", "order"])
def test_mesh2mesh_restart(launcher_script, request, tmp_path, change):
    """An interpolated restart matches a run on the target discretisation."""
    coarse = _genmeshbox(tmp_path, 4, "coarse.nmsh")
    fine = _genmeshbox(tmp_path, 6, "fine.nmsh")

    if change == "mesh":
        target_mesh, target_order = fine, 5
    else:
        target_mesh, target_order = coarse, 7

    # Keep the restart path short enough for Neko's fixed-size log buffer.
    dirs = {}
    for key in ("source", "reference", "restart"):
        context = TemporaryDirectory(prefix="neko-m2m-", dir="/tmp")
        request.addfinalizer(context.cleanup)
        dirs[key] = Path(context.name)

    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    name = request.node.name

    # Source run on the coarse mesh, writing a checkpoint after 2 steps.
    case = _case(coarse, dirs["source"], 5)
    case["case"]["time"]["end_time"] = CHECKPOINT_STEPS * TIMESTEP
    case["case"]["output_checkpoints"] = True
    case["case"]["checkpoint_control"] = "tsteps"
    case["case"]["checkpoint_value"] = CHECKPOINT_STEPS
    _run(
        launcher_script,
        case,
        tmp_path / "source.case",
        log_dir / f"{name}_source.log",
    )

    # Uninterrupted run on the target discretisation.
    case = _case(target_mesh, dirs["reference"], target_order)
    _run(
        launcher_script,
        case,
        tmp_path / "reference.case",
        log_dir / f"{name}_reference.log",
    )

    # Restart on the target discretisation.
    case = _case(target_mesh, dirs["restart"], target_order)
    case["case"]["restart_file"] = str(dirs["source"] / "fluid00000.chkp")
    if change == "mesh":
        case["case"]["restart_mesh_file"] = str(coarse)
    restart_log = log_dir / f"{name}_restart.log"
    _run(launcher_script, case, tmp_path / "restart.case", restart_log)

    assert "restarts at first order" in restart_log.read_text(), (
        "The interpolated restart did not reset the time integration order"
    )

    reference = _read_fld(dirs["reference"] / "field0.f00000")
    restarted = _read_fld(dirs["restart"] / "field0.f00000")
    for key in ("u", "v", "w", "s"):
        difference = np.max(np.abs(restarted[key] - reference[key]))
        assert difference < MAX_DIFFERENCE, (
            f"{key} after the interpolated restart differs from the "
            f"reference run by {difference}"
        )
