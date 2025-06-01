from dataclasses import dataclass
import pytest
import os
from os.path import join
from testlib import get_neko, get_neko_dir, run_neko
import json5

neko_dir = get_neko_dir()
examples_dir = join(neko_dir, "examples")
neko = get_neko()

@dataclass
class NekoTestCase:
    case_file: str = None
    user_file: str = None
    mesh_file: str = None


examples = {
    "hemi": NekoTestCase(
        case_file=join(examples_dir, "hemi", "hemi.case"),
        mesh_file=join(examples_dir, "hemi", "hemi.nmsh")
        ),
    "cylinder": NekoTestCase(
        case_file=join(examples_dir, "cylinder", "cylinder.case"),
        mesh_file=join(examples_dir, "cylinder", "cylinder.nmsh")
        ),
}

def manipulate_case(example, case):
    """Change the end_time in the case file to be twice the timestep and set the
    mesh file.

    Parameters
    ----------
    example : str
        The name of the example to modify.
    case : dict
        The dictionary representation of the case file.
    
    """

    case_object = case["case"]
    timestep = case_object.get("timestep", case_object.get("max_timestep"))
    case_object["end_time"] = 2 * timestep
    case_object["mesh_file"] = examples[example].mesh_file


def test_hemi(launcher_script, request, log_file):

    # Number of ranks to launch on
    nprocs = 2

    test_name = request.node.name
    case_file = examples["hemi"].case_file
    with open(case_file, "r") as f:
        case = json5.load(f)

    manipulate_case("hemi", case)

    case_file = join("tests", "test_examples", test_name + ".case") 

    with open(case_file, "w") as f:
        json.dump(case, f, indent=4)

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"


def test_cylinder(launcher_script, request):

    # Get the path to the neko executable

    # Get the name of the test from the request object
    test_name = request.node.name

    # Number of ranks to launch on
    nprocs = 2

    case_file = examples["cylinder"].case_file
    with open(case_file, "r") as f:
        case = json5.load(f)

    manipulate_case("cylinder", case)

    case_file = join("tests", "test_examples", test_name + ".case") 

    with open(case_file, "w") as f:
        json.dump(case, f, indent=4)

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"