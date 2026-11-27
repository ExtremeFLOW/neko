"""Integration coverage for all case-file users of ``import_fields``.

Flow and scalar initial conditions, sponge baseflows, and non-normal boundary
data all expose field import through the case file.  Exercise each pathway on
the original mesh and through global interpolation.  Interpolation is checked
both with coordinates embedded in the data file and with an explicitly named
mesh field.

ALE base-shape import is the only other production caller.  It does not expose
interpolation and its same-mesh path is covered by ``test_ale.py``.
"""

import json
import shutil
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

COMPONENTS = ("fluid", "scalar", "non_normal", "sponge")
IMPORT_MODES = ("same_mesh", "embedded_mesh", "separate_mesh")
NPROCS = 2
INTERPOLATION_TOLERANCE = {"dp": 1.0e-8, "sp": 1.0e-4}
VALUE_TOLERANCE = {"dp": 1.0e-8, "sp": 2.0e-4}


def _solver(solver_type):
    """Return a small, robust linear-solver configuration."""
    return {
        "type": solver_type,
        "preconditioner": {"type": "jacobi"},
        "absolute_tolerance": 1.0e-7,
        "max_iterations": 100,
    }


def _field_options(field_file, mode, mesh_field):
    """Build the shared case-file field import options."""
    options = {
        "file_name": str(field_file),
        "interpolate": mode != "same_mesh",
    }
    if mode == "separate_mesh":
        options["mesh_file_name"] = str(mesh_field)
    if mode != "same_mesh":
        options["interpolation"] = {
            "tolerance": INTERPOLATION_TOLERANCE[conftest.RP],
            "padding": 0.05,
        }
    return options


def _fluid(initial_condition, boundary_conditions, source_terms=None):
    """Build the fluid section shared by source and consumer cases."""
    fluid = {
        "scheme": "pnpn",
        "rho": 1.0,
        "mu": 1.0,
        "freeze": True,
        "initial_condition": initial_condition,
        "velocity_solver": _solver("cg"),
        "pressure_solver": _solver("gmres"),
        "boundary_conditions": boundary_conditions,
        "output_control": "never",
    }
    if source_terms is not None:
        fluid["source_terms"] = source_terms
    return fluid


def _scalar(initial_condition):
    """Build a passive scalar that preserves its imported constant value."""
    return {
        "enabled": True,
        "name": "temperature",
        "advection": False,
        "lambda": 1.0,
        "cp": 1.0,
        "initial_condition": initial_condition,
        "solver": _solver("cg"),
        "boundary_conditions": [
            {
                "type": "neumann",
                "zone_indices": [1, 2, 3, 4, 5, 6],
                "flux": 0.0,
            }
        ],
    }


def _base_case(mesh, output_directory):
    """Return the inexpensive one-step case shared by all runs."""
    return {
        "version": 1,
        "case": {
            "mesh_file": str(mesh),
            "output_directory": str(output_directory),
            "output_at_end": False,
            "output_boundary": False,
            "output_checkpoints": False,
            "time": {"end_time": 1.0e-3, "timestep": 1.0e-3},
            "numerics": {
                "time_order": 1,
                "polynomial_order": 3,
                "dealias": False,
            },
        },
    }


def _probe(fields, output_file, coordinates):
    """Build a one-point probe used to verify imported values."""
    return {
        "type": "probes",
        "compute_control": "nsamples",
        "compute_value": 1,
        "append_output": False,
        "output_file": output_file,
        "fields": fields,
        "points": [{"type": "points", "coordinates": coordinates}],
    }


def _source_case(mesh, workdir):
    """Create velocity and scalar source data with known values."""
    config = _base_case(mesh, workdir)
    case = config["case"]
    case["output_at_end"] = True
    case["fluid"] = _fluid(
        {
            "type": "expression",
            "value": ["1 + x", "2 + y", "3 + z"],
        },
        [
            {
                "type": "velocity_value",
                "zone_indices": [1, 2, 3, 4, 5, 6],
                "value": [1.5, 2.5, 3.5],
            }
        ],
    )
    case["fluid"]["output_filename"] = "source"
    case["fluid"]["output_mesh_in_all_files"] = True
    case["scalar"] = _scalar({"type": "uniform", "value": 4.0})
    return config


def _consumer_case(component, mode, assets, run_dir):
    """Create one consumer case and its expected probe values."""
    mesh = assets["source_mesh"] if mode == "same_mesh" else assets["target_mesh"]
    options = _field_options(assets["field"], mode, assets["mesh_field"])
    config = _base_case(mesh, run_dir)
    case = config["case"]

    no_slip = {
        "type": "no_slip",
        "zone_indices": [1, 2, 3, 4, 5, 6],
    }
    fluid_ic = {"type": "uniform", "value": [0.0, 0.0, 0.0]}
    boundaries = [no_slip]
    source_terms = None

    if component == "fluid":
        fluid_ic = {"type": "field", **options}
        probe_fields = ["u", "v", "w"]
        probe_point = [0.5, 0.5, 0.5]
        expected = [1.5, 2.5, 3.5]
    elif component == "scalar":
        # The first scalar in fluid output occupies the temperature slot,
        # whose import index is zero.
        scalar_ic = {"type": "field", "target_index": 0, **options}
        case["scalar"] = _scalar(scalar_ic)
        probe_fields = ["temperature"]
        probe_point = [0.5, 0.5, 0.5]
        expected = [4.0]
    elif component == "non_normal":
        boundaries = [
            {"type": "normal_outflow", "zone_indices": [2], **options},
            {"type": "no_slip", "zone_indices": [1, 3, 4, 5, 6]},
        ]
        probe_fields = ["u", "v", "w"]
        probe_point = [1.0, 0.5, 0.5]
        expected = [0.0, 2.5, 3.5]
    elif component == "sponge":
        source_terms = [
            {
                "type": "sponge",
                "amplitudes": [0.0, 0.0, 0.0],
                "fringe_registry_name": "u",
                "baseflow": {"method": "field", **options},
            }
        ]
        probe_fields = ["sponge_bf_u", "sponge_bf_v", "sponge_bf_w"]
        probe_point = [0.5, 0.5, 0.5]
        expected = [1.5, 2.5, 3.5]
    else:  # pragma: no cover - protects additions to COMPONENTS
        raise ValueError(f"Unknown import component: {component}")

    case["fluid"] = _fluid(fluid_ic, boundaries, source_terms)
    if component == "non_normal":
        # A frozen fluid does not apply boundary conditions during the time
        # loop. Take one non-advecting step so the imported tangential values
        # are observable on the boundary probe.
        case["fluid"]["freeze"] = False
        case["fluid"]["advection"] = False
    probe_name = f"probe_{component}_{mode}.csv"
    case["simulation_components"] = [_probe(probe_fields, probe_name, probe_point)]
    return config, run_dir / probe_name, expected


def _resolve_mesh_checker(neko):
    """Locate mesh_checker next to Neko or on PATH."""
    checker = which_command("mesh_checker")
    if checker is not None:
        return Path(checker).resolve()
    candidate = Path(neko).resolve().with_name("mesh_checker")
    assert candidate.is_file(), "The mesh_checker executable could not be found"
    return candidate


def _generate_mesh(genmeshbox, mesh_checker, workdir, name, dimensions):
    """Generate and validate one non-periodic box mesh."""
    result = subprocess.run(
        [
            str(genmeshbox),
            "0",
            "1",
            "0",
            "1",
            "0",
            "1",
            *(str(value) for value in dimensions),
            ".false.",
            ".false.",
            ".false.",
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout

    generated = workdir / "box.nmsh"
    mesh = workdir / name
    generated.rename(mesh)
    result = subprocess.run(
        [str(mesh_checker), mesh.name],
        cwd=workdir,
        capture_output=True,
        text=True,
        errors="replace",
    )
    assert result.returncode == 0, result.stdout
    return mesh


@pytest.fixture(scope="module")
def import_assets(tmp_path_factory, request):
    """Generate meshes and a single source field shared by the matrix."""
    workdir = tmp_path_factory.mktemp("import_fields")
    neko = Path(get_neko()).resolve()
    genmeshbox = Path(get_genmeshbox()).resolve()
    mesh_checker = _resolve_mesh_checker(neko)
    launcher = Path(request.config.getoption("--launcher-script")).resolve()

    source_mesh = _generate_mesh(
        genmeshbox, mesh_checker, workdir, "source.nmsh", (2, 2, 2)
    )
    target_mesh = _generate_mesh(
        genmeshbox, mesh_checker, workdir, "target.nmsh", (3, 2, 2)
    )

    source_case = workdir / "source.case"
    source_case.write_text(
        json.dumps(_source_case(source_mesh, workdir), indent=2) + "\n",
        encoding="utf-8",
    )
    source_log = workdir / "source.log"
    result = run_neko(
        str(launcher),
        configure_nprocs(NPROCS),
        str(source_case),
        str(neko),
        str(source_log),
    )
    assert result.returncode == 0, source_log.read_text(encoding="utf-8")

    field = workdir / "source0.f00000"
    assert field.is_file(), "The source run did not write source0.f00000"
    mesh_field = workdir / "source_mesh0.f00000"
    shutil.copyfile(field, mesh_field)

    return {
        "workdir": workdir,
        "neko": neko,
        "launcher": launcher,
        "source_mesh": source_mesh,
        "target_mesh": target_mesh,
        "field": field,
        "mesh_field": mesh_field,
    }


@pytest.mark.parametrize("component", COMPONENTS)
@pytest.mark.parametrize("mode", IMPORT_MODES)
def test_case_file_import_fields(component, mode, import_assets):
    """Import known data through every case-file pathway and verify values."""
    run_dir = import_assets["workdir"] / f"{component}_{mode}"
    run_dir.mkdir()
    config, probe_file, expected = _consumer_case(
        component, mode, import_assets, run_dir
    )
    case_file = run_dir / "consumer.case"
    case_file.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    log_file = run_dir / "consumer.log"

    result = run_neko(
        str(import_assets["launcher"]),
        configure_nprocs(NPROCS),
        str(case_file),
        str(import_assets["neko"]),
        str(log_file),
    )
    log = log_file.read_text(encoding="utf-8")
    assert result.returncode == 0, log
    assert "Import fields" in log
    assert f"Interpolation : {'T' if mode != 'same_mesh' else 'F'}" in log
    if mode == "separate_mesh":
        assert f"Mesh file     : {import_assets['mesh_field']}" in log

    values = [
        float(value) for value in probe_file.read_text().splitlines()[-1].split(",")
    ]
    assert values[1:] == pytest.approx(expected, abs=VALUE_TOLERANCE[conftest.RP])
