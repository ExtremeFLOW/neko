from os.path import join
from testlib import get_neko, run_neko, configure_nprocs, get_makeneko
import json5
import json
import subprocess


def test_user_stats(launcher_script, request, log_file, tmp_path):
    """A demo test. Notice that the parameters here are passed in via
    fixtures defined in conftest.py."""

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

    # Number of ranks to launch on
    max_nprocs = 1

    nprocs = configure_nprocs(max_nprocs)

    case_file = join("tests", "test_user_stats", "test_user_stats.case")

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

