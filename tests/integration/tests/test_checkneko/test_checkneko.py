"""Integration tests for checkneko, which checks case files without Neko.

Most tests edit a copy of the cylinder example. Its mesh has the labeled zones
1, 2, 5, 6 and 7, and its fluid boundary conditions use each of them once.
"""
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import json5
import pytest

from testlib import get_checkneko

REPO_ROOT = Path(__file__).resolve().parents[4]
EXAMPLES = REPO_ROOT / "examples"
CASE_FILES = (
    sorted(EXAMPLES.rglob("*.case"))
    + sorted(EXAMPLES.rglob("*.json"))
    + sorted((REPO_ROOT / "tests").rglob("*.case"))
)


def run(*args, cwd=None):
    """Run checkneko with the given arguments."""
    return subprocess.run(
        [str(Path(get_checkneko()).resolve()), *map(str, args)],
        cwd=cwd,
        capture_output=True,
        text=True,
        errors="replace",
    )


def problems(result):
    """The error and warning lines of a report."""
    return re.findall(r"^\S.*?: (?:error|warning): .*$", result.stdout, re.M)


def write_and_check(workdir, data, name="test.case"):
    """Write a case file into a directory and check it there."""
    (workdir / name).write_text(json.dumps(data, indent=4))
    return run(name, cwd=workdir)


@pytest.fixture
def cylinder(tmp_path):
    """The case of the cylinder example, with its mesh copied to tmp_path."""
    shutil.copy(EXAMPLES / "cylinder" / "cylinder.nmsh", tmp_path)
    with (EXAMPLES / "cylinder" / "cylinder.case").open() as handle:
        return json5.load(handle)


@pytest.mark.parametrize(
    "case_file",
    CASE_FILES,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_examples_pass_the_schema(case_file):
    """Every shipped case is valid, and every schema keyword is supported.

    Other checks may fail, since some examples generate their mesh first.
    """
    result = run(case_file, cwd=case_file.parent)
    assert re.search(r"Schema \.+ ok", result.stdout), \
        result.stdout + result.stderr
    assert "checkneko: warning" not in result.stderr, result.stderr


@pytest.mark.parametrize(
    "property_path",
    [
        ("case", "load_balance"),
        ("case", "joblimit"),
        ("case", "output_format"),
        ("case", "fluid", "variable_material_properties"),
    ],
)
def test_rejects_unsupported_keys(tmp_path, property_path):
    """Keys Neko ignores are rejected, as by validate_case_schema.py."""
    with (EXAMPLES / "tgv" / "tgv.case").open() as handle:
        data = json5.load(handle)
    target = data
    for key in property_path[:-1]:
        target = target[key]
    target[property_path[-1]] = False

    result = write_and_check(tmp_path, data)
    assert result.returncode == 1
    assert any("is not supported here" in p for p in problems(result))


def zone_table(result):
    """The rows of the zone table of a report, as lists of cells."""
    lines = result.stdout.splitlines()
    start = next(i for i, line in enumerate(lines)
                 if re.match(r"^\s*Zone\s+Facets\b", line))
    rows = []
    for line in lines[start:]:
        if not line.strip():
            break
        rows.append(re.split(r"\s{2,}", line.strip()))
    return rows


def test_cylinder_is_clean(tmp_path, cylinder):
    """The example itself has no problems."""
    result = write_and_check(tmp_path, cylinder)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "No problems found" in result.stdout
    assert "zones 1 2 5 6 7" in result.stdout
    assert zone_table(result) == [
        ["Zone", "Facets", "Fluid"],
        ["1", "72", "velocity_value"],
        ["2", "72", "outflow+dong"],
        ["5", "120", "normal_outflow"],
        ["6", "120", "normal_outflow"],
        ["7", "96", "no_slip"],
    ]


def set_zones(index, zones):
    """Edit that sets the zones of a fluid boundary condition."""
    def edit(case):
        case["fluid"]["boundary_conditions"][index]["zone_indices"] = zones
    return edit


def drop_bc(index):
    """Edit that removes a fluid boundary condition."""
    def edit(case):
        del case["fluid"]["boundary_conditions"][index]
    return edit


def drop_bcs(case):
    """Edit that removes all fluid boundary conditions."""
    del case["fluid"]["boundary_conditions"]


@pytest.mark.parametrize(
    "edit, message",
    [
        (set_zones(0, [21]), "the zone index 21 is out of range"),
        (set_zones(0, [1, 9]), "there is no zone 9 in the mesh"),
        (set_zones(1, [2, 7]),
         "zone 7 already has a fluid boundary condition"),
        (drop_bc(4), "the mesh zone 7 has no fluid boundary condition"),
        (drop_bcs, "there is no boundary_conditions entry"),
    ],
)
def test_zone_errors(tmp_path, cylinder, edit, message):
    """Zones of boundary conditions are checked against the mesh."""
    edit(cylinder["case"])
    result = write_and_check(tmp_path, cylinder)
    assert result.returncode == 1
    assert any(message in p for p in problems(result)), result.stdout


def test_compressible_zones(tmp_path):
    """The compressible solver takes several entries for one zone and leaves
    zones without one free, here the outflow, zone 2."""
    shutil.copy(EXAMPLES / "cyl_boundary_layer" / "cyl.nmsh", tmp_path)
    with (EXAMPLES / "euler_2d_cylinder" / "euler_2d_cylinder.case").open() \
            as handle:
        data = json5.load(handle)
    result = write_and_check(tmp_path, data)
    assert result.returncode == 0, result.stdout + result.stderr
    assert problems(result) == []
    rows = zone_table(result)
    assert ["1", "153", "velocity_value, density_value, pressure_value"] \
        in rows
    assert ["2", "153", "-"] in rows

    data["case"]["fluid"]["boundary_conditions"][3]["zone_indices"] = [9]
    result = write_and_check(tmp_path, data)
    assert result.returncode == 0
    assert any("warning:" in p and "there is no zone 9 in the mesh" in p
               for p in problems(result))


def test_scalars(tmp_path, cylinder):
    """Each scalar has a column in the zone table, named as by Neko, and the
    names must differ. A zone that is not in the mesh has a row too."""
    case = cylinder["case"]
    bcs = case["fluid"]["boundary_conditions"]
    bcs[4]["name"] = "cylinder"
    bcs[4]["zone_indices"] = "walls"
    case["constants"] = [{"name": "walls", "value": [7.0, 9.0]}]
    scalar_bcs = [
        {"type": "dirichlet", "value": 1.0, "zone_indices": [1]},
        {"type": "neumann", "flux": 0.0, "zone_indices": [2, 5, 6]},
        {"type": "dirichlet", "value": 0.0, "zone_indices": [7]},
    ]
    scalar = {
        "enabled": True, "Pe": 100.0,
        "initial_condition": {"type": "uniform", "value": 0.0},
        "solver": {"type": "cg", "preconditioner": {"type": "jacobi"},
                   "absolute_tolerance": 1e-6, "max_iterations": 100},
    }
    case["scalars"] = [
        dict(scalar, name="temperature", boundary_conditions=scalar_bcs),
        dict(scalar, boundary_conditions=scalar_bcs[:2]),
        dict(scalar, name="s_2", boundary_conditions=scalar_bcs),
    ]

    result = write_and_check(tmp_path, cylinder)
    assert result.returncode == 1
    rows = zone_table(result)
    assert rows[0] == ["Zone", "Facets", "Fluid", "Scalar temperature",
                       "Scalar s_2", "Scalar s_2"]
    assert ["7", "96", "no_slip (cylinder)", "dirichlet", "-",
            "dirichlet"] in rows
    assert ["9", "-", "no_slip (cylinder)", "-", "-", "-"] in rows
    found = problems(result)
    for message in [
        "there is no zone 9 in the mesh",
        "the mesh zone 7 has no scalar boundary condition",
        'the scalar name "s_2" is also the name of case.scalars[2]',
    ]:
        assert any(message in p for p in found), result.stdout


def test_constants(tmp_path, cylinder):
    """Values that name a constant need one of the kind Neko reads."""
    case = cylinder["case"]
    case["constants"] = [
        {"name": "inlet", "value": [1.0]},
        {"name": "U", "value": 1},
        {"name": "nu", "value": 0.01},
        {"name": "vel", "value": [1.0, 0.0]},
    ]
    case["fluid"]["boundary_conditions"][0]["zone_indices"] = "inlet"
    case["fluid"]["boundary_conditions"][0]["value"] = "vel"
    case["fluid"]["Re"] = "U"
    case["time"]["timestep"] = "dtt"
    case["numerics"]["polynomial_order"] = "nu"

    result = write_and_check(tmp_path, cylinder)
    assert result.returncode == 1
    found = problems(result)
    for message in [
        'the constant "dtt" is not defined in case.constants',
        'the constant "nu" must be an integer here',
        'the constant "U" must be a real number here',
        'the constant "vel" must have 3 entries here',
    ]:
        assert any(message in p for p in found), result.stdout
    # The zones of the inlet come from the constant
    assert not any("zone 1" in p for p in found), result.stdout


def test_point_zone_reference(tmp_path, cylinder):
    """A point zone a value refers to must be defined."""
    case = cylinder["case"]
    case["point_zones"] = [{
        "name": "box1", "geometry": "box", "x_bounds": [0.0, 1.0],
        "y_bounds": [0.0, 1.0], "z_bounds": [0.0, 1.0],
    }]
    case["fluid"]["initial_condition"] = {
        "type": "point_zone", "zone_name": "bx1",
        "base_value": [0.0, 0.0, 0.0], "zone_value": [1.0, 0.0, 0.0],
    }
    result = write_and_check(tmp_path, cylinder)
    assert result.returncode == 1
    assert any('the point zone "bx1" is not defined' in p
               and 'did you mean "box1"' in p for p in problems(result))


def test_missing_input_file(tmp_path, cylinder):
    """Files the case reads must exist."""
    cylinder["case"]["restart_file"] = "fluid00042.chkp"
    result = write_and_check(tmp_path, cylinder)
    assert result.returncode == 1
    assert any('the file "fluid00042.chkp" does not exist' in p
               for p in problems(result))


def test_truncated_mesh(tmp_path, cylinder):
    """A mesh shorter than its header says is reported. The file of the
    example is longer than that, Neko's writer leaves old bytes at its end."""
    data = (tmp_path / "cylinder.nmsh").read_bytes()
    (tmp_path / "cylinder.nmsh").write_bytes(data[:200000])
    result = write_and_check(tmp_path, cylinder)
    assert result.returncode == 1
    assert any("the file is truncated" in p for p in problems(result))


def test_misspelled_values_and_keys(tmp_path, cylinder):
    """Misspellings get a suggestion: an error where only known keys are
    allowed, a warning where any key is."""
    case = cylinder["case"]
    case["fluid"]["boundary_conditions"][4]["type"] = "no_slp"
    case["output_checkpoint"] = True
    result = write_and_check(tmp_path, cylinder)
    assert result.returncode == 1
    found = problems(result)
    assert any("error:" in p and 'got "no_slp" (did you mean "no_slip"?)' in p
               for p in found), result.stdout
    assert any("warning:" in p and 'did you mean "output_checkpoints"?' in p
               for p in found), result.stdout


def test_line_numbers_and_comments(tmp_path):
    """Comments are accepted as by Neko, and problems have the right line."""
    (tmp_path / "test.case").write_text("""\
// A comment
{
    "version": 1.0, # another one
    "case": {
        "mesh_file": "missing.nmsh", ! and a third
        "output_at_end": true,
        "output_at_end": false,
        "time": {"end_time": 1.0, "timestep": 1},
    }
}
""")
    result = run("test.case", cwd=tmp_path)
    assert result.returncode == 1
    found = problems(result)
    assert any(p.startswith("test.case:4: error: case: missing the required "
                            'key "numerics"') for p in found), result.stdout
    assert any(p.startswith("test.case:5: error: case.mesh_file:")
               for p in found), result.stdout
    assert any(p.startswith("test.case:7: warning: case.output_at_end:")
               for p in found), result.stdout
    assert any(p.startswith("test.case:8: error: case.time.timestep: "
                            "expected a real number") for p in found), \
        result.stdout


def test_parse_error(tmp_path):
    """Invalid JSON is reported with its line."""
    (tmp_path / "test.case").write_text(
        '{\n  "version": 1.0,\n  "case": {\n    "mesh_file": "x.nmsh",,\n'
        '  }\n}\n')
    result = run("test.case", cwd=tmp_path)
    assert result.returncode == 1
    assert "test.case:4: error: invalid JSON" in result.stdout


def test_schema_lookup(tmp_path, cylinder):
    """The schemas are found in share/neko/schemas next to the bin directory,
    also through a symbolic link on PATH, as from a Spack view."""
    (tmp_path / "test.case").write_text(json.dumps(cylinder, indent=4))
    prefix = tmp_path / "prefix"
    (prefix / "bin").mkdir(parents=True)
    exe = prefix / "bin" / "checkneko"
    shutil.copy(Path(get_checkneko()).resolve(), exe)

    def check(command, env=None):
        return subprocess.run(command + ["test.case"], cwd=tmp_path, env=env,
                              capture_output=True, text=True,
                              errors="replace")

    result = check([str(exe)])
    assert result.returncode == 1
    assert "cannot find the schemas" in result.stderr

    (prefix / "share" / "neko").mkdir(parents=True)
    (prefix / "share" / "neko" / "schemas").symlink_to(
        REPO_ROOT / "doc" / "schemas")
    result = check([str(exe)])
    assert result.returncode == 0, result.stdout + result.stderr

    view = tmp_path / "view" / "bin"
    view.mkdir(parents=True)
    (view / "checkneko").symlink_to(exe)
    env = dict(os.environ)
    env["PATH"] = str(view) + os.pathsep + env.get("PATH", "")
    result = check(["checkneko"], env=env)
    assert result.returncode == 0, result.stdout + result.stderr


def test_command_line(tmp_path, cylinder):
    """Options in the style of getopt(3)."""
    result = run()
    assert result.returncode == 1
    assert result.stderr.startswith("usage: checkneko")

    result = run("-h")
    assert result.returncode == 0
    assert result.stdout.startswith("usage: checkneko")

    result = run("-x", "test.case")
    assert result.returncode == 1
    assert "illegal option -- x" in result.stderr

    result = run("-d")
    assert result.returncode == 1
    assert "option requires an argument -- d" in result.stderr

    # Relative file names are relative to the run directory, by default the
    # current one
    (tmp_path / "test.case").write_text(json.dumps(cylinder, indent=4))
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    result = run(tmp_path / "test.case", cwd=elsewhere)
    assert result.returncode == 1
    # Without the mesh, the facets of the zones are not known
    assert ["7", "?", "no_slip"] in zone_table(result)
    result = run("-d", tmp_path, tmp_path / "test.case", cwd=elsewhere)
    assert result.returncode == 0, result.stdout + result.stderr
    result = run(f"-d{tmp_path}", "-s", REPO_ROOT / "doc" / "schemas",
                 tmp_path / "test.case", cwd=elsewhere)
    assert result.returncode == 0, result.stdout + result.stderr
