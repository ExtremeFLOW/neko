from os.path import join
from testlib import get_neko, run_neko, configure_nprocs, get_makeneko
import subprocess
import numpy as np
import json
import os


def read_rdcode(fld_path):
    """Return the 10-character rdcode from the NEKTON fld header."""
    with open(fld_path, "rb") as handle:
        header = handle.read(132)
    header_text = header.decode("ascii", errors="ignore")
    # rdcode is written by format '#std ... 1x,10a' which ends at column 93.
    # Use fixed positions to avoid trimming away trailing spaces in rdcode.
    if len(header_text) < 93:
        return header_text[-10:]
    return header_text[83:93]


def expected_rdcode_field(skip_pressure, skip_temperature):
    """Expected rdcode for writing a single field_t with mesh."""
    codes = ["X"]
    if skip_pressure and skip_temperature:
        codes.extend(["S", "0", "1"])
    elif skip_pressure:
        codes.append("T")
    else:
        codes.append("P")
    return ("".join(codes)).ljust(10)


def expected_rdcode_list(n_items, skip_pressure, skip_temperature, skip_velocity):
    """Expected rdcode for writing a field_list_t with mesh."""
    codes = ["X"]
    idx = 1

    write_pressure = not skip_pressure and idx <= n_items
    if write_pressure:
        idx += 1

    write_velocity = not skip_velocity and idx + 2 <= n_items
    if write_velocity:
        idx += 3

    write_temperature = not skip_temperature and idx <= n_items
    if write_temperature:
        idx += 1

    if write_velocity:
        codes.append("U")

    if write_pressure:
        codes.append("P")

    if write_temperature:
        codes.append("T")

    if idx <= n_items:
        n_scalars = n_items - idx + 1
        codes.append("S")
        codes.append(str(n_scalars // 10))
        codes.append(str(n_scalars % 10))

    return ("".join(codes)).ljust(10)


def test_fld_file_skip_params(launcher_script, request, log_file):
    """
    Verify rdcode generation for a single field with skip flag combinations.
    """
    test_name = request.node.name
    neko = "./neko"
    makeneko = get_makeneko()

    result = subprocess.run(
        [makeneko, join("tests", "test_fld_file", "test_fld_file.f90")],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    assert (
        result.returncode == 0
    ), f"makeneko process failed with exit code {result.returncode}"

    test_dir = join("tests", "test_fld_file")
    result = subprocess.run(
        ["genmeshbox", "0", "1", "0", "1", "0", "1", "3", "3", "3",
         ".true.", ".true.", ".true."],
        cwd=test_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    assert (
        result.returncode == 0
    ), f"genmeshbox process failed with exit code {result.returncode}"

    nprocs = configure_nprocs(1)

    base_case_path = join(test_dir, "test_fld_file.json")
    with open(base_case_path, "r", encoding="utf-8") as handle:
        case_data = json.load(handle)

    case_data["case"]["mesh_file"] = join(
        "tests", "test_fld_file", "box.nmsh"
    )

    fld_path = join(test_dir, "test0.f00000")

    cases = [
        (False, False, False),
        (True, False, False),
        (False, True, False),
        (False, False, True),
        (True, True, False),
        (True, False, True),
        (False, True, True),
        (True, True, True),
    ]

    for idx, (skip_p, skip_t, skip_v) in enumerate(cases):
        case_data["skip_pressure"] = skip_p
        case_data["skip_temperature"] = skip_t
        case_data["skip_velocity"] = skip_v
        case_data["is_list"] = False
        case_data["n_items"] = 1

        case_file = join(test_dir, f"test_fld_file_{idx}.json")
        with open(case_file, "w", encoding="utf-8") as handle:
            json.dump(case_data, handle, indent=2)
            handle.write("\n")

        if os.path.exists(fld_path):
            os.remove(fld_path)

        result = run_neko(
            launcher_script, nprocs, case_file, neko, log_file
        )
        assert (
            result.returncode == 0
        ), f"neko process failed with exit code {result.returncode}"

        assert os.path.exists(fld_path), f"{fld_path} was not created"
        rdcode = read_rdcode(fld_path)
        expected_rdcode = expected_rdcode_field(skip_p, skip_t)
        assert (
            rdcode == expected_rdcode
        ), f"rdcode mismatch for skip flags {skip_p, skip_t, skip_v}: {rdcode}"


def test_fld_file_skip_params_list_3(launcher_script, request, log_file):
    """
    Verify rdcode generation for a 3-item field_list_t with skip flags.
    """
    neko = "./neko"
    makeneko = get_makeneko()

    result = subprocess.run(
        [makeneko, join("tests", "test_fld_file", "test_fld_file.f90")],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    assert (
        result.returncode == 0
    ), f"makeneko process failed with exit code {result.returncode}"

    test_dir = join("tests", "test_fld_file")
    result = subprocess.run(
        ["genmeshbox", "0", "1", "0", "1", "0", "1", "3", "3", "3",
         ".true.", ".true.", ".true."],
        cwd=test_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    assert (
        result.returncode == 0
    ), f"genmeshbox process failed with exit code {result.returncode}"

    nprocs = configure_nprocs(1)

    base_case_path = join(test_dir, "test_fld_file.json")
    with open(base_case_path, "r", encoding="utf-8") as handle:
        case_data = json.load(handle)

    case_data["case"]["mesh_file"] = join(
        "tests", "test_fld_file", "box.nmsh"
    )

    fld_path = join(test_dir, "test0.f00000")

    cases = [
        (False, False, False),
        (True, False, False),
        (False, True, False),
        (False, False, True),
        (True, True, False),
        (True, False, True),
        (False, True, True),
        (True, True, True),
    ]

    for idx, (skip_p, skip_t, skip_v) in enumerate(cases):
        case_data["skip_pressure"] = skip_p
        case_data["skip_temperature"] = skip_t
        case_data["skip_velocity"] = skip_v
        case_data["is_list"] = True
        case_data["n_items"] = 3

        case_file = join(test_dir, f"test_fld_file_list3_{idx}.json")
        with open(case_file, "w", encoding="utf-8") as handle:
            json.dump(case_data, handle, indent=2)
            handle.write("\n")

        if os.path.exists(fld_path):
            os.remove(fld_path)

        result = run_neko(
            launcher_script, nprocs, case_file, neko, log_file
        )
        assert (
            result.returncode == 0
        ), f"neko process failed with exit code {result.returncode}"

        assert os.path.exists(fld_path), f"{fld_path} was not created"
        rdcode = read_rdcode(fld_path)
        expected_rdcode = expected_rdcode_list(3, skip_p, skip_t, skip_v)
        assert (
            rdcode == expected_rdcode
        ), f"rdcode mismatch for skip flags {skip_p, skip_t, skip_v}: {rdcode}"
