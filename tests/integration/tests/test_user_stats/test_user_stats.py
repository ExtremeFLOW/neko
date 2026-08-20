from os.path import join
from testlib import get_neko, run_neko, configure_nprocs, get_makeneko, get_genmeshbox
import subprocess
import json
import numpy as np


def test_user_stats(launcher_script, request, log_file, tmp_path):
    """Runs a simple case where 1 scalar is assigned random values from [0,1].
    We test whether the mean after averaging across xy is close to 0.5."""

    # Gets the nameof the test, i.e. test_demo here. `request` can be use for
    # other things like this.
    test_name = request.node.name

    # Get the path to the neko executable
    neko = "./neko"
    makeneko = get_makeneko()

    result = subprocess.run(
        [makeneko, join("tests", "test_user_stats", "test_user_stats.f90")],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True)

    assert (
        result.returncode == 0
    ), f"makeneko process failed with exit code {result.returncode}"

    genmeshbox = get_genmeshbox()
    result = subprocess.run(
        [genmeshbox, "0", "1", "0", "1", "0", "1", "3", "3", "3",
         ".true.", ".true.", ".true."],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True)
    assert (
        result.returncode == 0
    ), f"genmeshbox process failed with exit code {result.returncode}"

    # Number of ranks to launch on
    max_nprocs = 1

    nprocs = configure_nprocs(max_nprocs)

    case_template = join("tests", "test_user_stats", "test_user_stats.case")
    with open(case_template, "r", encoding="utf-8") as handle:
        case_data = json.load(handle)
    case_data["case"]["output_directory"] = str(tmp_path)

    case_file = tmp_path / "test_user_stats.case"
    with open(case_file, "w", encoding="utf-8") as handle:
        json.dump(case_data, handle, indent=2)
        handle.write("\n")

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

    stats_file = tmp_path / "user_stats0.csv"
    assert stats_file.exists(), (
        f"user_stats output was not written to {tmp_path}"
    )

    csv = np.genfromtxt(stats_file, delimiter=",")
    mean = np.mean(csv[12:25, 2])
    error = (mean - 0.5) / 0.5

    assert (
         error < 1e-4
    ), f"Error exceeded tolerance: {error}"
