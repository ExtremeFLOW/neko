import subprocess
import pytest
import os


def get_neko():
    """
    Returns the path to the neko executable  via the
    environmental variable NEKO_BIN, or just returns "neko" 
    """

    src_dir = os.getenv("NEKO_BIN")
    if src_dir:
        return os.path.join(src_dir, "neko")
    else:
        return "neko"


def test_cylinder():
    """A very good description of the test."""

    # Get the path to the neko executable
    neko = get_neko()

    # Number of ranks to launch on
    nprocs = 2

    log_file = "cylinder.log"

    # Either specify a case, or load it here into a json, manipulate
    # and save a new case file.
    case_file = "cylinder.case"

    # Run Neko
    with open("logs/" + log_file, "w") as log:
        result = subprocess.run(
            ["mpirun", "-n", str(nprocs), neko, f"case_templates/{case_file}"],
            stdout=log,
            stderr=subprocess.STDOUT,
        )

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

    assert (False), "Intentional failure"
