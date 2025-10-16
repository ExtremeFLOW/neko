"""
Tests using a 3D cylinder case.
Checks
- Deformations.
- Residuals.
- Restarts for the fluid.
"""

from os.path import join
import os
from testlib import get_neko, run_neko, configure_nprocs, get_neko_dir, parse_log
import json5
import json
import numpy as np
from numpy.testing import assert_allclose
from conftest import RP


def test_cylinder(launcher_script, request, tmp_path):
    """
    First part of the test, which just runs the case for a bit.

    Note, we handle log names manually here to have _part1 and part2.
    The reason we do everything in a single test is that we want to compare
    logs from 2 runs, one with a restart and one without.

    For single precision, the tests is much more relaxed, we reduce the
    tolerance a lot and only check pressure residuals.
    """

    # Make sure logs directory exists
    os.makedirs("logs", exist_ok=True)

    # Set the precision for the test
    eps = 1e-15 if RP == "dp" else 1e-2

    test_name = request.node.name + "_part1"
    log_file = os.path.join("logs", f"{test_name}.log")

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
    ref_data = np.genfromtxt(join(test_dir, "reflog1_" + RP + ".csv"),
         delimiter=",", names=True)

    # Parse the log file from the test run
    parsed_log = join(test_dir, "log1.csv")
    parse_log(log_file, parsed_log)

    # Load the parsed data
    parsed_data_part1 = np.genfromtxt(parsed_log, delimiter=",", names=True)

    # Compare
    columns = parsed_data_part1.dtype.names # Get column names read from  the CSV

    for i in columns:
        if i == "total_step_time" or (RP == "sp" and "velocity" in i):
            continue
        assert_allclose(parsed_data_part1[i], ref_data[i], rtol=eps,
                        err_msg=f"Column '{i}' does not match reference data.")

    #
    # Part 2, run from chekpoint file
    #

    test_name = request.node.name + "_part2"
    log_file = os.path.join("logs", f"{test_name}.log")

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
    ref_data = np.genfromtxt(join(test_dir, "reflog2_" + RP + ".csv"),
        delimiter=",", names=True)

    # Parse the log file from the test run
    parse_log(log_file, join(test_dir, "log2.csv"))

    # Load the parsed data
    parsed_data_part2 = np.genfromtxt(join(test_dir, "log2.csv"), delimiter=",",
        names=True)

    for i in columns:
        if i == "total_step_time" or (RP == "sp" and "velocity" in i):
            # Skip total_step_time as it may vary
            continue
        assert_allclose(parsed_data_part2[i], ref_data[i], rtol=eps,
                        err_msg=f"Column '{i}' does not match reference data.")

    #
    # Part 3, compare overlapping time steps from the 2 runs
    #

    parsed_data_part1 = parsed_data_part1[parsed_data_part1["time"] > 5e-2]
    for i in columns:
        if i in ["total_step_time", "step"] or (RP == "sp" and "velocity" in i):
            continue
        assert_allclose(parsed_data_part2[i], parsed_data_part1[i], rtol=eps,
                        err_msg=f"Column '{i}' does not match reference data.")

    # Harsher tolerance on the first overlapped iteration even for SP.
    assert (
        (parsed_data_part1["pressure_start_residual"][0] -
         parsed_data_part2["pressure_start_residual"][0]) / parsed_data_part1["pressure_start_residual"][0] < 1e-6
    )

