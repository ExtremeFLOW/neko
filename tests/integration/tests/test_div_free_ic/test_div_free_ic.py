"""Integration coverage for the divergence-free initial-condition projection.

``case.fluid.initial_condition.make_divergence_free`` projects the initial
velocity onto the divergence-free subspace before the first time step.  The
projection has a closed-form solution for the case set up here, so the test
compares against it rather than against previously recorded output.

On the unit cube, with the velocity prescribed at ``x = 0`` and on the walls at
``y = 0, 1``, the pressure prescribed at ``x = 1`` and periodicity in ``z``, the
initial condition

    u = 1 + 0.3 sin(2 pi x),  v = 0.3 sin(2 pi y),  w = 0

is the gradient part of a Helmholtz--Leray decomposition plus

    u = 1 - 0.3 cos(2 pi y) sinh(2 pi x) / cosh(2 pi),  v = w = 0,

which is the field the projection must produce.  Note that it leaves the
prescribed inflow at ``x = 0`` untouched while changing the interior, so the
comparison also pins down the boundary treatment of the projection.
"""

import json
import math
import subprocess
from pathlib import Path

import pytest

import conftest
from testlib import (
    configure_nprocs,
    get_genmeshbox,
    get_neko,
    run_neko,
    which_command,
)

NPROCS = 2

# Number of elements per direction and the polynomial order of the test mesh.
# At least three elements are needed in the periodic direction, a mesh with two
# does not get a valid periodic connectivity.
ELEMENTS = (4, 4, 3)
POLYNOMIAL_ORDER = 7

# Probe locations, all on the centre line y = z = 0.5.
PROBE_X = (0.0, 0.5, 1.0)

# Tolerance of the comparison against the closed-form projection.
VALUE_TOLERANCE = {"dp": 1.0e-5, "sp": 1.0e-3}

# The divergence of the initial condition is sqrt(0.36 pi^2) and has to come
# down by at least this factor.
DIVERGENCE_REDUCTION = 1.0e4


def _analytic_u(x):
    """The x velocity of the projected field on the centre line y = 0.5."""
    return 1.0 + 0.3 * math.sinh(2.0 * math.pi * x) / math.cosh(2.0 * math.pi)


def _solver(solver_type, preconditioner):
    """Return the linear-solver configuration used for both runs."""
    return {
        "type": solver_type,
        "preconditioner": {"type": preconditioner},
        "absolute_tolerance": 1.0e-10,
        "max_iterations": 500,
    }


def _case(mesh, output_directory, project):
    """Build the case file, with the projection either on or off."""
    initial_condition = {
        "type": "expression",
        "value": ["1 + 0.3*sin(2*pi*x)", "0.3*sin(2*pi*y)", "0"],
    }
    if project:
        initial_condition["make_divergence_free"] = True

    coordinates = []
    for x in PROBE_X:
        coordinates += [x, 0.5, 0.5]

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
                "polynomial_order": POLYNOMIAL_ORDER,
                "dealias": False,
            },
            "fluid": {
                "scheme": "pnpn",
                "Re": 100.0,
                # Freeze the fluid so the probes see the initial condition, and
                # not a field that has been advanced by a time step.
                "freeze": True,
                "initial_condition": initial_condition,
                "boundary_conditions": [
                    {
                        "type": "velocity_value",
                        "zone_indices": [1],
                        "value": [1.0, 0.0, 0.0],
                    },
                    {"type": "outflow", "zone_indices": [2]},
                    {"type": "no_slip", "zone_indices": [3, 4]},
                ],
                "velocity_solver": _solver("cg", "jacobi"),
                "pressure_solver": _solver("gmres", "hsmg"),
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
                    "points": [{"type": "points", "coordinates": coordinates}],
                }
            ],
        },
    }


def _resolve_mesh_checker(neko):
    """Locate mesh_checker next to Neko or on PATH."""
    checker = which_command("mesh_checker")
    if checker is not None:
        return Path(checker).resolve()
    candidate = Path(neko).resolve().with_name("mesh_checker")
    assert candidate.is_file(), "The mesh_checker executable could not be found"
    return candidate


def _generate_mesh(genmeshbox, mesh_checker, workdir):
    """Generate a unit cube, periodic in z only, and validate it."""
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

    mesh = workdir / "box.nmsh"
    result = subprocess.run(
        [str(mesh_checker), mesh.name],
        cwd=workdir,
        capture_output=True,
        text=True,
        errors="replace",
    )
    assert result.returncode == 0, result.stdout
    return mesh


def _read_probes(probe_file):
    """Return the last sampled u, v, w triplet for each probe point."""
    rows = [
        line.strip()
        for line in probe_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    # The file opens with a header and one line of coordinates per point, then
    # carries one line of values per point and sample.
    samples = rows[1 + len(PROBE_X):]
    assert len(samples) >= len(PROBE_X), rows
    values = []
    for row in samples[-len(PROBE_X):]:
        # Each row is time, u, v, w.
        values.append([float(entry) for entry in row.split(",")][1:])
    return values


def _run(assets, project, name):
    """Run one case and return the probed values and the log."""
    run_dir = assets["workdir"] / name
    run_dir.mkdir()
    case_file = run_dir / f"{name}.case"
    case_file.write_text(
        json.dumps(_case(assets["mesh"], run_dir, project), indent=2) + "\n",
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
    return _read_probes(run_dir / "probe.csv"), log


@pytest.fixture(scope="module")
def div_free_assets(tmp_path_factory, request):
    """Generate the mesh shared by both runs."""
    workdir = tmp_path_factory.mktemp("div_free_ic")
    neko = Path(get_neko()).resolve()
    genmeshbox = Path(get_genmeshbox()).resolve()
    launcher = Path(request.config.getoption("--launcher-script")).resolve()
    mesh = _generate_mesh(genmeshbox, _resolve_mesh_checker(neko), workdir)
    return {
        "workdir": workdir,
        "neko": neko,
        "launcher": launcher,
        "mesh": mesh,
    }


def test_divergence_free_initial_condition(div_free_assets):
    """The projected initial condition matches the closed-form solution."""
    tolerance = VALUE_TOLERANCE[conftest.RP]

    unprojected, _ = _run(div_free_assets, False, "plain")
    projected, log = _run(div_free_assets, True, "projected")

    # Without the keyword the initial condition is used as given. All three
    # probes sit at sin(2 pi x) = 0, so the x velocity is exactly one there.
    for x, values in zip(PROBE_X, unprojected):
        assert values[0] == pytest.approx(1.0, abs=tolerance), f"x = {x}"
        assert values[1] == pytest.approx(0.0, abs=tolerance), f"x = {x}"
        assert values[2] == pytest.approx(0.0, abs=tolerance), f"x = {x}"

    # With the keyword the field is the divergence-free one, which differs from
    # the initial condition everywhere except on the inflow boundary.
    for x, values in zip(PROBE_X, projected):
        assert values[0] == pytest.approx(_analytic_u(x), abs=tolerance), \
            f"x = {x}"
        assert values[1] == pytest.approx(0.0, abs=tolerance), f"x = {x}"
        assert values[2] == pytest.approx(0.0, abs=tolerance), f"x = {x}"

    # The projection has to be visible, i.e. the interior really did change.
    assert abs(projected[2][0] - 1.0) > 0.1

    # And the reported divergence has to come down by orders of magnitude.
    before = _log_value(log, "div(u) L2, before :")
    after = _log_value(log, "div(u) L2, after  :")
    assert after * DIVERGENCE_REDUCTION < before, f"{before} -> {after}"


def _log_value(log, marker):
    """Extract the number reported after a marker in the Neko log."""
    for line in log.splitlines():
        if marker in line:
            return float(line.split(marker)[1].strip())
    raise AssertionError(f"'{marker}' is missing from the log")
