import subprocess
import pytest
import os
from testlib import get_neko


def test_cylinder(launcher_script, request):
    """A very good description of the test."""

    # Get the path to the neko executable
    neko = get_neko()

    # Get the name of the test from the request object
    test_name = request.node.name


    # Number of ranks to launch on
    nprocs = 2

    log_file = "logs/" + test_name + ".log"

    # Either specify a case, or load it here into a json, manipulate
    # and save a new case file.
    case_file = "cylinder.case"

    # Run Neko
    with open(log_file, "w") as log:
        result = subprocess.run(
            [launcher_script, str(nprocs), case_file, neko],
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True
        )

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

