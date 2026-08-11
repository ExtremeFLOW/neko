"""Test scalar solver residuals and restart continuity."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from numpy.testing import assert_allclose

from conftest import BACKEND, RP
from testlib import configure_nprocs, get_neko, get_neko_dir, parse_log, run_neko


SCALAR_TOLERANCE = 1e-8
if BACKEND == "cpu" and RP == "dp":
    REFERENCE_RTOL = 1e-12
    REFERENCE_ATOL = 1e-16
elif BACKEND == "cpu":
    REFERENCE_RTOL = 1e-6
    REFERENCE_ATOL = 1e-12
# Looser for the device because there is no clean reference data
elif RP == "dp":
    REFERENCE_RTOL = 1e-7
    REFERENCE_ATOL = 1e-12
else:
    REFERENCE_RTOL = 1e-4
    REFERENCE_ATOL = 1e-12
if RP == "dp":
    OVERLAP_RTOL = 1e-11
    OVERLAP_ATOL = 1e-13
else:
    OVERLAP_RTOL = 1e-6
    OVERLAP_ATOL = 1e-12


def _load_log(log_file, parsed_file):
    """Parse a Neko log and return its tabular data."""
    expected_precision = "double precision" if RP == "dp" else "single precision"
    assert expected_precision in log_file.read_text(encoding="utf-8"), (
        f"Neko executable precision does not match --real_precision={RP}"
    )
    parse_log(str(log_file), str(parsed_file))
    return np.genfromtxt(parsed_file, delimiter=",", names=True)


def _compare_with_reference(actual, reference):
    """Compare deterministic scalar convergence data with a reference run."""
    for column in (
        "time",
        "dt",
        "tracer_iters",
        "tracer_start_residual",
        "tracer_final_residual",
    ):
        assert_allclose(
            actual[column],
            reference[column],
            rtol=REFERENCE_RTOL,
            atol=REFERENCE_ATOL,
            err_msg=f"Column '{column}' does not match reference data.",
        )

    assert np.all(actual["tracer_final_residual"] <= SCALAR_TOLERANCE)


def test_scalar_restart(launcher_script, request, tmp_path):
    """Check scalar residuals and continuity across a checkpoint restart."""
    neko = get_neko()
    neko_dir = Path(get_neko_dir()).resolve()
    test_dir = (
        neko_dir / "tests" / "integration" / "tests" / "test_scalar_restart"
    )
    nprocs = configure_nprocs(1)
    # Keep the restart path short enough for Neko's fixed-size log buffer.
    checkpoint_dir_context = TemporaryDirectory(
        prefix="neko-scalar-", dir="/tmp"
    )
    request.addfinalizer(checkpoint_dir_context.cleanup)
    checkpoint_dir = Path(checkpoint_dir_context.name)
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)

    with (test_dir / "scalar_restart.case").open(encoding="utf-8") as stream:
        case = json.load(stream)

    case["case"]["mesh_file"] = str(
        neko_dir / "tests" / "integration" / "meshes" / "small_test_cyl.nmsh"
    )
    case["case"]["output_directory"] = str(checkpoint_dir)

    part1_case = tmp_path / "scalar_restart_part1.case"
    with part1_case.open("w", encoding="utf-8") as stream:
        json.dump(case, stream, indent=2)

    part1_log = log_dir / f"{request.node.name}_part1.log"
    result = run_neko(
        launcher_script, nprocs, str(part1_case), neko, str(part1_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )

    part1_data = _load_log(part1_log, tmp_path / "part1.csv")
    part1_reference = np.genfromtxt(
        test_dir / f"reference_part1_{BACKEND}_{RP}.csv", delimiter=",", names=True
    )
    _compare_with_reference(part1_data, part1_reference)

    checkpoints = sorted(checkpoint_dir.glob("scalar_restart*.chkp"))
    assert len(checkpoints) == 3, (
        "Expected initial, t=0.05, and t=0.10 checkpoints, "
        f"found {checkpoints}"
    )

    case["case"]["restart_file"] = str(checkpoints[1])
    case["case"]["output_checkpoints"] = False
    del case["case"]["checkpoint_control"]
    del case["case"]["checkpoint_value"]

    part2_case = tmp_path / "scalar_restart_part2.case"
    with part2_case.open("w", encoding="utf-8") as stream:
        json.dump(case, stream, indent=2)

    part2_log = log_dir / f"{request.node.name}_part2.log"
    result = run_neko(
        launcher_script, nprocs, str(part2_case), neko, str(part2_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )

    part2_data = _load_log(part2_log, tmp_path / "part2.csv")
    part2_reference = np.genfromtxt(
        test_dir / f"reference_part2_{BACKEND}_{RP}.csv", delimiter=",", names=True
    )
    _compare_with_reference(part2_data, part2_reference)

    overlap = part1_data[part1_data["time"] > 0.05]
    for column in (
        "time",
        "dt",
        "tracer_iters",
        "tracer_start_residual",
        "tracer_final_residual",
    ):
        assert_allclose(
            part2_data[column],
            overlap[column],
            rtol=OVERLAP_RTOL,
            atol=OVERLAP_ATOL,
            err_msg=f"Restart changed overlapping column '{column}'.",
        )

    assert np.any(part1_data["tracer_start_residual"] > 0.0)
