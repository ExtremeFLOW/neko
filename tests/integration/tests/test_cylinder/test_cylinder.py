"""
Tests using a 3D cylinder case.
Checks
- Deformations.
- Residuals.
- Restarts for the fluid.
"""

from os.path import join
from testlib import get_neko, run_neko, configure_nprocs, get_neko_dir, parse_log
import json5
import json
import numpy as np
from numpy.testing import assert_allclose


def test_cylinder_part1(launcher_script, request, log_file, tmp_path):
    """
    First part of the test, which just runs the case for a bit.
    """

    # Gets the nameof the test, i.e. test_demo here. `request` can be use for
    # other things like this.
    test_name = request.node.name

    # Get the path to the neko executable
    neko = get_neko()
    neko_dir = get_neko_dir()
    test_dir = join(neko_dir, "tests", "integration", "tests", "test_cylinder")

    max_nprocs = 1

    nprocs = configure_nprocs(max_nprocs)

    # Run Neko
    result = run_neko(launcher_script, nprocs, join(test_dir,
        "cylinder_part1.case"), neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

    # Reference parsed log data
    ref_data = np.genfromtxt(join(test_dir, "reflog1_dp.csv"), delimiter=",",
        names=True)

    # Parse the log file from the test run
    parse_log(log_file, join(test_dir, "log1_dp.csv"))

    # Load the parsed data
    parsed_data = np.genfromtxt(join(test_dir, "log1_dp.csv"), delimiter=",",
        names=True)

    # Compare
    columns = parsed_data.dtype.names # Get column names read from  the CSV

    for i in columns:
        if i == "total_step_time":
            # Skip total_step_time as it may vary
            continue
        assert_allclose(parsed_data[i], ref_data[i], rtol=1e-15,
                        err_msg=f"Column '{i}' does not match reference data.")

def test_cylinder_part2(launcher_script, request, log_file, tmp_path):
    """
    Second part of the test, which test the restart functionality.
    """

    # Gets the nameof the test, i.e. test_demo here. `request` can be use for
    # other things like this.
    test_name = request.node.name

    # Get the path to the neko executable
    neko = get_neko()
    neko_dir = get_neko_dir()
    test_dir = join(neko_dir, "tests", "integration", "tests", "test_cylinder")

    max_nprocs = 1

    nprocs = configure_nprocs(max_nprocs)

    # We start with the p1 file and moidfy it
    case_file_p1 = join(test_dir, "cylinder_part1.case")

    with open(case_file_p1, "r") as f:
        case = json.load(f)

    case["case"]["restart_file"] = "fluid00001.chkp"

    case_file = join("tests", "test_demo", test_name + ".case")
    with open(case_file, "w") as f:
        json.dump(case, f, indent=4)

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

    # Reference parsed log data
    ref_data = np.genfromtxt(join(test_dir, "reflog2_dp.csv"), delimiter=",",
        names=True)

    # Parse the log file from the test run
    parse_log(log_file, join(test_dir, "log2_dp.csv"))

    # Load the parsed data
    parsed_data = np.genfromtxt(join(test_dir, "log2_dp.csv"), delimiter=",",
        names=True)

    # Compare
    columns = parsed_data.dtype.names # Get column names read from  the CSV

    for i in columns:
        if i == "total_step_time":
            # Skip total_step_time as it may vary
            continue
        assert_allclose(parsed_data[i], ref_data[i], rtol=1e-15,
                        err_msg=f"Column '{i}' does not match reference data.")
