from os.path import join
from testlib import get_neko, run_neko, configure_nprocs, get_neko_dir
import json5
import json


def test_demo(launcher_script, request, log_file, tmp_path):
    """A demo test. Notice that the parameters here are passed in via
    fixtures defined in conftest.py."""

    # Gets the nameof the test, i.e. test_demo here. `request` can be use for
    # other things like this.
    test_name = request.node.name

    # Get the path to the neko executable
    neko = get_neko()
    neko_dir = get_neko_dir()

    # Number of ranks to launch on
    max_nprocs = 2

    # nprocs = min(max_nprocs, MAX_NPROCS)
    nprocs = configure_nprocs(max_nprocs)

    # Either specify a case, or load it here into a json, manipulate
    # and save a new case file.
    case_file = join(neko_dir, "examples", "cylinder", "cylinder.case")

    # Read the case file. We use json5 to allow comments in the case file.
    with open(case_file, "r") as f:
        case = json5.load(f)

    # Manipulate the case file
    case_object = case["case"]
    case_object["mesh_file"] = join("meshes", "small_test_cyl.nmsh")
    # Fix the end time to be twice the timestep
    time_object = case_object["time"]
    timestep = time_object.get("timestep", time_object.get("max_timestep"))
    time_object["end_time"] = 2 * timestep
    # Use the tmp_path fixture to set the output directory to a temporary path,
    # so that one does not pollute the working directory with output files.
    case_object["output_directory"] = str(tmp_path)

    case_file = join("tests", "test_demo", test_name + ".case")
    with open(case_file, "w") as f:
        json.dump(case, f, indent=4)


    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

