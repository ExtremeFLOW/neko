from copy import deepcopy
import json
from pathlib import Path
from testlib import run_neko, configure_nprocs, get_makeneko
import subprocess


def _build_and_mesh():
    """Build user file and generate mesh; reuse across tests."""
    neko = "./neko"
    makeneko = get_makeneko()

    build = subprocess.run(
        [makeneko, "tests/test_scalars_ic/test_scalars_ic.f90"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    assert build.returncode == 0, (
        f"makeneko failed with exit code {build.returncode}\n{build.stdout}"
    )

    mesh = subprocess.run(
        [
            "genmeshbox",
            "0",
            "1",
            "0",
            "1",
            "0",
            "1",
            "3",
            "3",
            "3",
            ".true.",
            ".true.",
            ".true.",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    assert mesh.returncode == 0, (
        f"genmeshbox failed with exit code {mesh.returncode}\n{mesh.stdout}"
    )
    return neko


def test_scalars_ic_single(launcher_script, log_file, tmp_path):
    """Run scalar initial conditions for each IC variant in ics.json (single scalar)."""

    neko = _build_and_mesh()
    nprocs = configure_nprocs(1)

    base_case = json.loads(Path("tests/test_scalars_ic/test_scalars_ic.json").read_text())
    ic_variants = json.loads(Path("tests/test_scalars_ic/ics.json").read_text())

    for name, ic_cfg in ic_variants.items():
        case = deepcopy(base_case)
        case["case"]["scalar"]["initial_condition"] = ic_cfg

        case_file = Path("tests/test_scalars_ic") / f"test_scalars_ic_single_{name}.case"
        case_file.write_text(json.dumps(case, indent=2))

        result = run_neko(launcher_script, nprocs, str(case_file), neko, log_file)
        assert result.returncode == 0, (
            f"neko failed for single-scalar IC '{name}' with exit code {result.returncode}"
        )


def test_scalars_ic_multiple(launcher_script, log_file):
    """Run scalar initial conditions for two scalars, varying temperature ICs."""

    neko = _build_and_mesh()
    nprocs = configure_nprocs(1)

    base_case = json.loads(Path("tests/test_scalars_ic/test_scalars_ic.json").read_text())
    scalars_array = json.loads(Path("tests/test_scalars_ic/scalars.json").read_text())
    ic_variants = json.loads(Path("tests/test_scalars_ic/ics.json").read_text())

    # Replace single scalar with scalars array
    case = deepcopy(base_case)
    # Remove legacy single-scalar block
    case["case"].pop("scalar", None)
    case["case"]["scalars"] = scalars_array["scalars"]

    for name, ic_cfg in ic_variants.items():
        case_variant = deepcopy(case)

        # Adjust temperature IC values for specific variants
        ic_variant = deepcopy(ic_cfg)
        if name == "uniform":
            ic_variant["value"] = 2.0
        elif name == "point_zone":
            ic_variant["base_value"] = 2.0
            ic_variant["zone_value"] = 2.0
        elif name == "field":
            ic_variant["file_name"] = "tests/test_scalars_ic/ics_multi.f00000"

        # Update only the temperature scalar's IC
        for scalar in case_variant["case"]["scalars"]:
            if scalar.get("name", "") == "temperature":
                scalar["initial_condition"] = ic_variant

        case_file = Path("tests/test_scalars_ic") / f"test_scalars_ic_multi_{name}.case"
        case_file.write_text(json.dumps(case_variant, indent=2))

        result = run_neko(launcher_script, nprocs, str(case_file), neko, log_file)
        assert result.returncode == 0, (
            f"neko failed for multi-scalar IC '{name}' with exit code {result.returncode}"
        )
