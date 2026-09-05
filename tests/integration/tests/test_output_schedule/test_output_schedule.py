"""Test which files an output schedule produces, and that a restart
continues it without losing or repeating a file.

The case writes both the fluid fields and a checkpoint every 0.05 time
units, and runs from 0 to 0.1 with a time step of 0.01, with
`output_at_end` on. Note that accumulating that time step ten times lands a
few ulps short of the end time, which used to buy the run an eleventh step
past `end_time` and, with it, a duplicate of every output.
"""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

from testlib import configure_nprocs, get_neko, get_neko_dir, run_neko


def _write_case(path, case, output_dir, mesh, restart_file=None):
    """Write the case file, pointing it at the working directory."""
    case = json.loads(json.dumps(case))
    case["case"]["mesh_file"] = str(mesh)
    case["case"]["output_directory"] = str(output_dir)
    if restart_file is not None:
        case["case"]["restart_file"] = str(restart_file)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(case, stream, indent=2)
    return path


def _files(directory, pattern):
    """The names of the files matching a pattern, in order."""
    return sorted(path.name for path in directory.glob(pattern))


def test_output_schedule(launcher_script, request, tmp_path):
    """Check the files written by a run and by a restart of it."""
    neko = get_neko()
    neko_dir = Path(get_neko_dir()).resolve()
    test_dir = (
        neko_dir / "tests" / "integration" / "tests" / "test_output_schedule"
    )
    mesh = neko_dir / "tests" / "integration" / "meshes" / "small_test_cyl.nmsh"
    nprocs = configure_nprocs(1)

    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)

    # Keep the restart path short enough for Neko's fixed-size log buffer.
    work_context = TemporaryDirectory(prefix="neko-out-", dir="/tmp")
    request.addfinalizer(work_context.cleanup)
    work_dir = Path(work_context.name)

    with (test_dir / "output_schedule.case").open(encoding="utf-8") as stream:
        case = json.load(stream)

    #
    # Part 1, the whole interval in one run
    #
    part1_dir = work_dir / "part1"
    part1_dir.mkdir()
    part1_case = _write_case(
        tmp_path / "part1.case", case, part1_dir, mesh
    )
    part1_log = log_dir / f"{request.node.name}_part1.log"
    result = run_neko(
        launcher_script, nprocs, str(part1_case), neko, str(part1_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )

    # The fields are written at t = 0, 0.05 and 0.1. The write forced at the
    # end of the run does not add a fourth file, since t = 0.1 is written by
    # the last time step already.
    assert _files(part1_dir, "field0.f*") == [
        "field0.f00000",
        "field0.f00001",
        "field0.f00002",
    ], f"Unexpected field files: {_files(part1_dir, 'field0.f*')}"

    # The checkpoints are written at t = 0.05 and 0.1. The initial condition
    # is not checkpointed, and the end of the run is not checkpointed twice.
    assert _files(part1_dir, "chkp*.chkp") == [
        "chkp00000.chkp",
        "chkp00001.chkp",
    ], f"Unexpected checkpoints: {_files(part1_dir, 'chkp*.chkp')}"

    #
    # Part 2, the second half of the same interval, restarted from the
    # checkpoint written at t = 0.05
    #
    part2_dir = work_dir / "part2"
    part2_dir.mkdir()
    part2_case = _write_case(
        tmp_path / "part2.case",
        case,
        part2_dir,
        mesh,
        restart_file=part1_dir / "chkp00000.chkp",
    )
    part2_log = log_dir / f"{request.node.name}_part2.log"
    result = run_neko(
        launcher_script, nprocs, str(part2_case), neko, str(part2_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )

    # Only what the first run had not written yet, numbered so that the two
    # runs together produce the files of the uninterrupted run: nothing is
    # written at the restart time, and nothing scheduled is skipped.
    assert _files(part2_dir, "field0.f*") == ["field0.f00002"], (
        f"Unexpected field files after the restart: "
        f"{_files(part2_dir, 'field0.f*')}"
    )
    assert _files(part2_dir, "chkp*.chkp") == ["chkp00001.chkp"], (
        f"Unexpected checkpoints after the restart: "
        f"{_files(part2_dir, 'chkp*.chkp')}"
    )

    #
    # Part 3, the same run with the initial state switched off
    #
    part3_dir = work_dir / "part3"
    part3_dir.mkdir()
    case_no_start = json.loads(json.dumps(case))
    case_no_start["case"]["output_at_start"] = False
    part3_case = _write_case(
        tmp_path / "part3.case", case_no_start, part3_dir, mesh
    )
    part3_log = log_dir / f"{request.node.name}_part3.log"
    result = run_neko(
        launcher_script, nprocs, str(part3_case), neko, str(part3_log)
    )
    assert result.returncode == 0, (
        f"neko process failed with exit code {result.returncode}"
    )

    # The fields are then written at t = 0.05 and 0.1 only.
    assert _files(part3_dir, "field0.f*") == [
        "field0.f00000",
        "field0.f00001",
    ], (
        f"Unexpected field files without the initial state: "
        f"{_files(part3_dir, 'field0.f*')}"
    )
