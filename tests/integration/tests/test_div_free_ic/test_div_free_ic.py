"""Integration coverage for the divergence-free initial-condition projection.

``case.fluid.initial_condition.make_divergence_free`` imposes the velocity
boundary conditions on the initial condition and projects it onto the
divergence-free subspace before the first time step. The projection has a
closed-form solution for the case set up here, so the test compares against it
rather than against previously recorded output.

On the unit cube, periodic in y and z, with the velocity prescribed at
``x = 0`` and the pressure at ``x = 1``, the initial condition

    u = 1 + 0.3 sin(2 pi x),  v = 0.3 sin(2 pi y),  w = 0

is the gradient part of a Helmholtz--Leray decomposition plus

    u = 1 - 0.3 cos(2 pi y) sinh(2 pi x) / cosh(2 pi),
    v = 0.3 sin(2 pi y) cosh(2 pi x) / cosh(2 pi),  w = 0,

which is the field the projection must produce; on the probed centre line
y = 0.5 its v vanishes. It leaves the prescribed inflow at ``x = 0``
untouched while changing the interior, so the comparison also pins down the
boundary treatment of the projection.
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
# At least three elements are needed in a periodic direction, a mesh with two
# does not get a valid periodic connectivity.
ELEMENTS = (4, 4, 3)
POLYNOMIAL_ORDER = 7

# Tolerance of the comparison against the closed-form projection, and the
# residual reduction the projection is solved to for it.
VALUE_TOLERANCE = {"dp": 1.0e-5, "sp": 1.0e-3}
SOLVE_TOLERANCE = {"dp": 1.0e-9, "sp": 1.0e-6}

# Required reduction of the divergence norm. With an inflow the prescribed
# tangential velocity there disagrees with the divergence-free interior, and
# the correction is masked on the boundary as in the time loop, so a one-point
# layer of divergence remains at the inflow and bounds the reduction. Without
# boundaries the reduction is only limited by the solve.
DIVERGENCE_REDUCTION = {"analytic": 10.0, "periodic": 1.0e4}

CENTRE_LINE = [(x, 0.5, 0.5) for x in (0.0, 0.5, 1.0)]
INFLOW = (0.0, 0.5, 0.5)
WALL = (0.5, 0.0, 0.5)


def _analytic_u(x):
    """The x velocity of the projected field on the centre line y = 0.5."""
    return 1.0 + 0.3 * math.sinh(2.0 * math.pi * x) / math.cosh(2.0 * math.pi)


def _solver(solver_type, preconditioner):
    """Return the linear-solver configuration used for all runs."""
    return {
        "type": solver_type,
        "preconditioner": {"type": preconditioner},
        "absolute_tolerance": 1.0e-10,
        "max_iterations": 500,
    }


def _case(mesh, output_directory, project, kind, points):
    """Build the case file for one of three initial conditions.

    ``analytic`` and ``uniform`` run on a mesh periodic in y and z with an
    inflow prescribing exactly the values of the initial condition, so that
    imposing the boundary conditions changes nothing and the closed form
    applies. ``periodic`` is a fully periodic box, a pure Neumann problem,
    whose projected field is the uniform (1, 0, 0). ``violating`` runs in a
    duct with walls, with an initial condition satisfying neither the inflow
    nor the wall condition.
    """
    if kind == "periodic":
        initial_condition = {
            "type": "expression",
            "value": ["1 + 0.3*sin(2*pi*x)", "0.3*sin(2*pi*y)",
                      "0.3*sin(2*pi*z)"],
        }
        boundary_conditions = None
    elif kind == "analytic":
        initial_condition = {
            "type": "expression",
            "value": ["1 + 0.3*sin(2*pi*x)", "0.3*sin(2*pi*y)", "0"],
        }
        boundary_conditions = [
            {"type": "expression_velocity", "zone_indices": [1],
             "value": ["1", "0.3*sin(2*pi*y)", "0"]},
            {"type": "outflow", "zone_indices": [2]},
        ]
    elif kind == "uniform":
        initial_condition = {"type": "uniform", "value": [1.0, 0.0, 0.0]}
        boundary_conditions = [
            {"type": "velocity_value", "zone_indices": [1],
             "value": [1.0, 0.0, 0.0]},
            {"type": "outflow", "zone_indices": [2]},
        ]
    elif kind == "violating":
        initial_condition = {"type": "uniform", "value": [0.7, 0.3, 0.0]}
        boundary_conditions = [
            {"type": "velocity_value", "zone_indices": [1],
             "value": [1.0, 0.0, 0.0]},
            {"type": "outflow", "zone_indices": [2]},
            {"type": "no_slip", "zone_indices": [3, 4]},
        ]
    else:  # pragma: no cover
        raise ValueError(kind)

    if project:
        initial_condition["make_divergence_free"] = True
        initial_condition["divergence_free_tolerance"] = \
            SOLVE_TOLERANCE[conftest.RP]
        initial_condition["divergence_free_max_iterations"] = 200

    case = {
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
                # Frozen, so the probes see the projected initial condition
                # and not a field advanced by a time step. A frozen fluid does
                # not impose its boundary conditions in the time loop either.
                "freeze": True,
                "initial_condition": initial_condition,
                "velocity_solver": _solver("cg", "jacobi"),
                # hsmg stalls on the singular pure-Neumann problem in single
                # precision, phmg does not.
                "pressure_solver": _solver(
                    "gmres", "phmg" if kind == "periodic" else "hsmg"),
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
                    "points": [{"type": "points",
                                "coordinates": [c for p in points for c in p]}],
                }
            ],
        },
    }
    if boundary_conditions is not None:
        case["case"]["fluid"]["boundary_conditions"] = boundary_conditions
    return case


def _resolve_mesh_checker(neko):
    """Locate mesh_checker next to Neko or on PATH."""
    checker = which_command("mesh_checker")
    if checker is not None:
        return Path(checker).resolve()
    candidate = Path(neko).resolve().with_name("mesh_checker")
    assert candidate.is_file(), "The mesh_checker executable could not be found"
    return candidate


def _generate_mesh(genmeshbox, mesh_checker, workdir, name, periodic):
    """Generate a unit cube with the given periodicity flags and validate it."""
    result = subprocess.run(
        [
            str(genmeshbox),
            "0", "1", "0", "1", "0", "1",
            *(str(value) for value in ELEMENTS),
            *(".true." if flag else ".false." for flag in periodic),
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout

    mesh = workdir / name
    (workdir / "box.nmsh").rename(mesh)
    result = subprocess.run(
        [str(mesh_checker), mesh.name],
        cwd=workdir,
        capture_output=True,
        text=True,
        errors="replace",
    )
    assert result.returncode == 0, result.stdout
    return mesh


def _read_probes(probe_file, npoints):
    """Return the last sampled (u, v, w) keyed by probe coordinates.

    The file opens with a header and one line of coordinates per point, then
    one line of ``time, u, v, w`` per point and sample. The points are not
    necessarily in the order they were requested, so the values are matched to
    the coordinates by position.
    """
    rows = [
        line.strip()
        for line in probe_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    coords = [tuple(round(float(e), 6) for e in r.split(","))
              for r in rows[1:1 + npoints]]
    samples = rows[1 + npoints:]
    assert len(samples) >= npoints, rows
    values = {}
    for point, row in zip(coords, samples[-npoints:]):
        values[point] = [float(e) for e in row.split(",")][1:]
    return values


def _run(assets, project, kind, points, name):
    """Run one case and return the probed values and the log."""
    run_dir = assets["workdir"] / name
    run_dir.mkdir()
    mesh = {"violating": assets["duct"],
            "periodic": assets["box"]}.get(kind, assets["periodic"])
    case_file = run_dir / f"{name}.case"
    case_file.write_text(
        json.dumps(_case(mesh, run_dir, project, kind, points), indent=2)
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
    return _read_probes(run_dir / "probe.csv", len(points)), log


def _div_norms(log):
    """Extract the divergence norm before and after from the Neko log."""
    marker = "||div u||_L2      :"
    for line in log.splitlines():
        if marker in line:
            before, after = line.split(marker)[1].split("->")
            if after.strip().startswith("unchanged"):
                after = before
            return float(before), float(after)
    raise AssertionError(f"'{marker}' is missing from the log")


@pytest.fixture(scope="module")
def div_free_assets(tmp_path_factory, request):
    """Generate the two meshes shared by the tests."""
    workdir = tmp_path_factory.mktemp("div_free_ic")
    neko = Path(get_neko()).resolve()
    genmeshbox = Path(get_genmeshbox()).resolve()
    mesh_checker = _resolve_mesh_checker(neko)
    launcher = Path(request.config.getoption("--launcher-script")).resolve()
    return {
        "workdir": workdir,
        "neko": neko,
        "launcher": launcher,
        "periodic": _generate_mesh(genmeshbox, mesh_checker, workdir,
                                   "periodic.nmsh", (False, True, True)),
        "duct": _generate_mesh(genmeshbox, mesh_checker, workdir,
                               "duct.nmsh", (False, False, True)),
        "box": _generate_mesh(genmeshbox, mesh_checker, workdir,
                              "box.nmsh", (True, True, True)),
    }


def test_divergence_free_initial_condition(div_free_assets):
    """The projected initial condition matches the closed-form solution."""
    tolerance = VALUE_TOLERANCE[conftest.RP]

    plain, _ = _run(div_free_assets, False, "analytic", CENTRE_LINE, "plain")
    projected, log = _run(div_free_assets, True, "analytic", CENTRE_LINE,
                          "projected")

    # Without the keyword the initial condition is used as given. All three
    # probes sit at sin(2 pi x) = 0, so the x velocity is exactly one there.
    for point in CENTRE_LINE:
        u, v, w = plain[point]
        assert u == pytest.approx(1.0, abs=tolerance), point
        assert v == pytest.approx(0.0, abs=tolerance), point
        assert w == pytest.approx(0.0, abs=tolerance), point

    # With the keyword the field is the divergence-free one, which differs from
    # the initial condition everywhere except on the inflow boundary.
    for point in CENTRE_LINE:
        u, v, w = projected[point]
        assert u == pytest.approx(_analytic_u(point[0]), abs=tolerance), point
        assert v == pytest.approx(0.0, abs=tolerance), point
        assert w == pytest.approx(0.0, abs=tolerance), point

    # The projection has to be visible, i.e. the interior really did change.
    assert abs(projected[CENTRE_LINE[2]][0] - 1.0) > 0.1

    # And the reported divergence has to come down.
    before, after = _div_norms(log)
    assert after * DIVERGENCE_REDUCTION["analytic"] < before, \
        f"{before} -> {after}"


def test_divergence_free_periodic_box(div_free_assets):
    """In a fully periodic box the divergent part is the whole perturbation,
    and the projection returns the uniform mean flow."""
    tolerance = VALUE_TOLERANCE[conftest.RP]
    points = [(0.37, 0.61, 0.25), (0.8, 0.2, 0.7)]

    projected, log = _run(div_free_assets, True, "periodic", points,
                          "periodic")
    for point in points:
        u, v, w = projected[point]
        assert u == pytest.approx(1.0, abs=tolerance), point
        assert v == pytest.approx(0.0, abs=tolerance), point
        assert w == pytest.approx(0.0, abs=tolerance), point

    before, after = _div_norms(log)
    assert after * DIVERGENCE_REDUCTION["periodic"] < before, \
        f"{before} -> {after}"
    # The pure Neumann path reports the net boundary flux, zero here.
    assert "Net boundary flux" in log


def test_divergence_free_leaves_a_good_field_alone(div_free_assets):
    """Projecting an already divergence-free field must not damage it."""
    tolerance = VALUE_TOLERANCE[conftest.RP]

    projected, log = _run(div_free_assets, True, "uniform", CENTRE_LINE,
                          "uniform")

    for point in CENTRE_LINE:
        u, v, w = projected[point]
        assert u == pytest.approx(1.0, abs=tolerance), point
        assert v == pytest.approx(0.0, abs=tolerance), point
        assert w == pytest.approx(0.0, abs=tolerance), point

    # The correction has to stay at the level of the initial divergence, the
    # projection must not introduce any of its own.
    before, after = _div_norms(log)
    assert after < max(1.0e3 * before, tolerance), f"{before} -> {after}"


def test_divergence_free_imposes_the_boundary_conditions(div_free_assets):
    """The velocity boundary conditions are imposed before the projection."""
    tolerance = VALUE_TOLERANCE[conftest.RP]
    points = [INFLOW, WALL]

    # Frozen and unprojected, the initial condition is left as it is, so the
    # probes see its values rather than those of the boundary conditions.
    plain, _ = _run(div_free_assets, False, "violating", points,
                    "violating_plain")
    assert plain[INFLOW][0] == pytest.approx(0.7, abs=tolerance)
    assert plain[WALL][1] == pytest.approx(0.3, abs=tolerance)

    projected, _ = _run(div_free_assets, True, "violating", points,
                        "violating")
    # The correction is masked on the strong velocity boundaries, so the
    # imposed values are kept exactly: u = 1 at the inflow, v = 0 at the wall.
    assert projected[INFLOW][0] == pytest.approx(1.0, abs=tolerance)
    assert projected[WALL][1] == pytest.approx(0.0, abs=tolerance)
