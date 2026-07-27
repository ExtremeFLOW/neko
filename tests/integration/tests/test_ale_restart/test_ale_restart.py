"""Test ALE checkpoint restart across polynomial orders."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from testlib import configure_nprocs, get_neko, get_neko_dir, run_neko


def _hdf5_enabled(neko_dir):
    """Return whether the configured Neko build includes HDF5 support."""
    makefile = neko_dir / "Makefile"
    if not makefile.exists():
        return False
    return "-DHAVE_HDF5=1" in makefile.read_text(encoding="utf-8")


def test_ale_hdf5_restart_different_polynomial_order(
    launcher_script, request, tmp_path
):
    """Restart HDF5 ALE data at a higher polynomial order."""
    neko = get_neko()
    neko_dir = Path(get_neko_dir()).resolve()
    if not _hdf5_enabled(neko_dir):
        pytest.skip("Neko was configured without HDF5 support")

    nprocs = configure_nprocs(2)
    checkpoint_dir_context = TemporaryDirectory(
        prefix="neko-ale-hdf5-", dir="/tmp"
    )
    request.addfinalizer(checkpoint_dir_context.cleanup)
    checkpoint_dir = Path(checkpoint_dir_context.name)
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)

    source_case = (
        neko_dir
        / "examples"
        / "oscillating_cylinder"
        / "oscillating_cylinder.case"
    )
    with source_case.open(encoding="utf-8") as stream:
        case = json.load(stream)

    case["case"]["mesh_file"] = str(
        neko_dir / "tests" / "integration" / "meshes" / "small_test_cyl.nmsh"
    )
    case["case"]["output_directory"] = str(checkpoint_dir)
    case["case"]["output_at_end"] = False
    case["case"]["output_boundary"] = False
    case["case"]["checkpoint_format"] = "hdf5"
    case["case"]["checkpoint_filename"] = "ale_restart"
    case["case"]["checkpoint_value"] = 1
    case["case"]["time"]["end_time"] = 0.002
    case["case"]["numerics"]["polynomial_order"] = 3
    case["case"]["fluid"]["output_control"] = "never"
    case["case"]["fluid"]["ale"]["solver"]["output_base_shape"] = False

    part1_case = tmp_path / "ale_hdf5_restart_part1.case"
    with part1_case.open("w", encoding="utf-8") as stream:
        json.dump(case, stream, indent=2)

    part1_log = log_dir / f"{request.node.name}_part1.log"
    result = run_neko(
        launcher_script, nprocs, str(part1_case), neko, str(part1_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )

    checkpoints = sorted(checkpoint_dir.glob("ale_restart*.h5"))
    assert len(checkpoints) >= 3, (
        "Expected initial and per-step ALE checkpoints, "
        f"found {checkpoints}"
    )

    case["case"]["restart_file"] = str(checkpoints[-2])
    case["case"]["output_checkpoints"] = False
    del case["case"]["checkpoint_control"]
    del case["case"]["checkpoint_value"]
    case["case"]["numerics"]["polynomial_order"] = 4

    part2_case = tmp_path / "ale_hdf5_restart_part2.case"
    with part2_case.open("w", encoding="utf-8") as stream:
        json.dump(case, stream, indent=2)

    part2_log = log_dir / f"{request.node.name}_part2.log"
    result = run_neko(
        launcher_script, nprocs, str(part2_case), neko, str(part2_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )
