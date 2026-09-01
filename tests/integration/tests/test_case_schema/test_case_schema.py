import json
from pathlib import Path
import re
import subprocess
import sys

import json5
import pytest


REPO_ROOT = Path(__file__).resolve().parents[4]
VALIDATOR = REPO_ROOT / "contrib" / "validate_case_schema.py"
EXAMPLES_DIR = REPO_ROOT / "examples"
TESTS_DIR = REPO_ROOT / "tests"
CASE_FILES = (
    sorted(EXAMPLES_DIR.rglob("*.case"))
    + sorted(EXAMPLES_DIR.rglob("*.json"))
    + sorted(TESTS_DIR.rglob("*.case"))
)


@pytest.mark.parametrize(
    "case_file",
    CASE_FILES,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_example_case_schema(case_file):
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"schema validation failed for {case_file}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )


@pytest.mark.parametrize(
    "property_path",
    [
        ("case", "load_balance"),
        ("case", "joblimit"),
        ("case", "output_format"),
        ("case", "fluid", "variable_material_properties"),
    ],
)
def test_rejects_known_ignored_properties(tmp_path, property_path):
    with (EXAMPLES_DIR / "tgv" / "tgv.case").open(encoding="utf-8") as handle:
        data = json5.load(handle)

    target = data
    for key in property_path[:-1]:
        target = target[key]
    target[property_path[-1]] = False

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


def test_allows_user_extension_properties(tmp_path):
    with (EXAMPLES_DIR / "tgv" / "tgv.case").open(encoding="utf-8") as handle:
        data = json5.load(handle)
    data["case"]["fluid"]["user_parameter"] = 42

    case_file = tmp_path / "custom.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_rejects_real_encoding_for_integer_controller_value(tmp_path):
    with (EXAMPLES_DIR / "rayleigh_benard_cylinder" / "rayleigh.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)
    data["case"]["fluid"]["output_value"] = 10.0

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


def test_rejects_integer_encoding_for_real_parameter(tmp_path):
    with (EXAMPLES_DIR / "api" / "c" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)
    data["case"]["fluid"]["Re"] = 200

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


def test_rejects_integer_encoding_in_real_array(tmp_path):
    with (EXAMPLES_DIR / "api" / "c" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)
    boundary = data["case"]["fluid"]["boundary_conditions"][0]
    boundary["value"] = [1, 0, 0]

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


@pytest.mark.parametrize(
    "boundary_type",
    ["outflow", "normal_outflow", "outflow+user", "normal_outflow+user"],
)
@pytest.mark.parametrize("property_name", ["delta", "velocity_scale"])
def test_rejects_dong_parameters_for_non_dong_outflow(
    tmp_path, boundary_type, property_name
):
    with (EXAMPLES_DIR / "cylinder" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)

    boundary = data["case"]["fluid"]["boundary_conditions"][1]
    boundary["type"] = boundary_type
    boundary[property_name] = 0.1

    case_file = tmp_path / "invalid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0


@pytest.mark.parametrize("boundary_type", ["outflow+dong", "normal_outflow+dong"])
def test_allows_dong_outflow_parameters(tmp_path, boundary_type):
    with (EXAMPLES_DIR / "cylinder" / "cylinder.case").open(
        encoding="utf-8"
    ) as handle:
        data = json5.load(handle)

    boundary = data["case"]["fluid"]["boundary_conditions"][1]
    boundary["type"] = boundary_type
    if boundary_type == "normal_outflow+dong":
        boundary["value"] = [1.0, 0.0, 0.0]
    boundary["delta"] = 0.01
    boundary["velocity_scale"] = 1.0

    case_file = tmp_path / "valid.case"
    case_file.write_text(json.dumps(data), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def collect_const_values(value, property_name):
    values = set()
    if isinstance(value, dict):
        property_schema = value.get(property_name)
        if isinstance(property_schema, dict) and "const" in property_schema:
            values.add(property_schema["const"])
        for child in value.values():
            values.update(collect_const_values(child, property_name))
    elif isinstance(value, list):
        for child in value:
            values.update(collect_const_values(child, property_name))
    return values


def factory_case_values(path):
    source = path.read_text(encoding="utf-8")
    return set(re.findall(r"case \([\"']([^\"']+)[\"']\)", source))


def referenced_definition_consts(schema, union_name, property_name):
    definition = schema["$defs"][union_name]
    property_schema = definition.get("properties", {}).get(property_name, {})
    if "const" in property_schema:
        return {property_schema["const"]}

    values = set()
    for variant in definition["oneOf"]:
        definition_name = variant["$ref"].removeprefix("#/$defs/")
        values.update(
            referenced_definition_consts(schema, definition_name, property_name)
        )
    return values


def test_all_builtin_simcomps_have_schemas():
    runtime_types = factory_case_values(
        REPO_ROOT / "src/simulation_components/simulation_component_fctry.f90"
    )
    with (REPO_ROOT / "doc/schemas/simulation-components.schema.json").open(
        encoding="utf-8"
    ) as handle:
        schema = json.load(handle)
    refs = json.dumps(schema["$defs"]["simulationComponent"])
    schema_types = set(re.findall(r"urn:neko:schema:simcomps:([^#]+)#", refs))

    assert schema_types == runtime_types


@pytest.mark.parametrize(
    ("factory", "schema", "property_name", "extra_runtime_types"),
    [
        (
            "src/wall_models/wall_model_fctry.f90",
            "doc/schemas/fluid-pnpn-boundary-conditions.schema.json",
            "model",
            set(),
        ),
        (
            "src/source_terms/source_term_fctry.f90",
            "doc/schemas/source-terms.schema.json",
            "type",
            {"user"},
        ),
        (
            "src/mesh/point_zone_fctry.f90",
            "doc/schemas/point-zones.schema.json",
            "geometry",
            {"combine"},
        ),
    ],
)
def test_all_builtin_factory_types_have_schemas(
    factory, schema, property_name, extra_runtime_types
):
    runtime_types = factory_case_values(REPO_ROOT / factory) | extra_runtime_types
    with (REPO_ROOT / schema).open(encoding="utf-8") as handle:
        schema_data = json.load(handle)
    if schema == "doc/schemas/source-terms.schema.json":
        schema_types = referenced_definition_consts(
            schema_data, "fluidSourceTerm", property_name
        )
    else:
        schema_types = collect_const_values(schema_data, property_name)

    assert schema_types == runtime_types
