from os.path import join
from testlib import get_neko, run_neko, configure_nprocs


def test_demo(launcher_script, log_file):
    """A demo test. Notice that the parameters here are passed in via 
    fixtures defined in conftest.py."""

    # Get the path to the neko executable
    neko = get_neko()

    # Number of ranks to launch on
    max_nprocs = 2

    nprocs = configure_nprocs(max_nprocs)
    # Either specify a case, or load it here into a json, manipulate
    # and save a new case file.
    case_file = join("case_templates", "cylinder.case")

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

