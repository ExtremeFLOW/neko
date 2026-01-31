from copy import deepcopy
import json
from pathlib import Path
from testlib import get_neko, run_neko, configure_nprocs, get_makeneko
import subprocess


def test_scalars_ic(launcher_script, log_file, tmp_path):
    """Run scalar initial conditions for each IC variant in ics.json."""

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

    nprocs = configure_nprocs(1)

    base_case = json.loads(Path("tests/test_scalars_ic/test_scalars_ic.json").read_text())
    ic_variants = json.loads(Path("tests/test_scalars_ic/ics.json").read_text())

    for name, ic_cfg in ic_variants.items():
        case = deepcopy(base_case)
        case["case"]["scalar"]["initial_condition"] = ic_cfg

        case_file = tmp_path / f"test_scalars_ic_{name}.case"
        case_file.write_text(json.dumps(case, indent=2))

        result = run_neko(launcher_script, nprocs, str(case_file), neko, log_file)
        assert result.returncode == 0, (
            f"neko failed for IC '{name}' with exit code {result.returncode}"
        )
