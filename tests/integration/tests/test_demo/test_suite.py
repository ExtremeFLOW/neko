import subprocess
import pytest
import os
from os.path import join
from testlib import get_neko, run_neko


def test_demo(launcher_script, log_file):
    """A very good description of the test."""

    # Get the path to the neko executable
    neko = get_neko()

    # Number of ranks to launch on
    nprocs = 2

    # Either specify a case, or load it here into a json, manipulate
    # and save a new case file.
    case_file = join("case_templates", "cylinder.case")

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

