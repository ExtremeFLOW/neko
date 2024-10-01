import subprocess
import pytest


def test_cylinder():
    """A very good description of the test."""

    log_file = "cylinder.log"

    # Either specify a case, or load it here into a json, manipulate
    # and save a new case file.
    case_file = "cylinder.case"

    # Run Neko
    with open("logs/" + log_file, "w") as log:
        result = subprocess.run(
            ["neko", f"case_templates/{case_file}"],
            stdout=log,
            stderr=subprocess.STDOUT,
        )

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"
