from dataclasses import dataclass
import subprocess
import pytest
import os
from os.path import join
from testlib import get_neko, get_neko_dir, run_neko, configure_nprocs, get_makeneko
import json5
import json

neko_dir = get_neko_dir()
examples_dir = join(neko_dir, "examples")
turboneko = get_neko()

"""A simple structure to hold the test case information for Neko examples."""
@dataclass
class NekoTestCase:
    case_file: str = None
    user_file: str = None
    mesh_file: str = None

examples_dir = join(neko_dir, "examples")

examples = {
    "2d_cylinder": NekoTestCase(
        case_file=join(examples_dir, "2d_cylinder", "2d_cylinder.case"),
        user_file=join(examples_dir, "2d_cylinder", "2d_cylinder.f90"),
        mesh_file=join(examples_dir, "2d_cylinder", "2d_cylinder.nmsh")
    ),
    "TS_channel": NekoTestCase(
        case_file=join(examples_dir, "TS_channel", "TS_channel.case"),
        user_file=join(examples_dir, "TS_channel", "TS_channel.f90"),
        mesh_file=join(examples_dir, "TS_channel", "box.nmsh")
    ),
    "advecting_cone": NekoTestCase(
        case_file=join(examples_dir, "advecting_cone", "advecting_cone.case"),
        user_file=join(examples_dir, "advecting_cone", "advecting_cone.f90"),
        mesh_file=join(examples_dir, "advecting_cone", "box.nmsh")
    ),
    "cyl_boundary_layer": NekoTestCase(
        case_file=join(examples_dir, "cyl_boundary_layer", "cyl_bl_basic.case"),
        mesh_file=join(examples_dir, "cyl_boundary_layer", "cyl.nmsh"),
        user_file=join(examples_dir, "cyl_boundary_layer", "cyl_bl.f90")
    ),
    "cylinder": NekoTestCase(
        case_file=join(examples_dir, "cylinder", "cylinder.case"),
        mesh_file=join(examples_dir, "cylinder", "cylinder.nmsh")
    ),
    "euler_1d_sod": NekoTestCase(
        case_file=join(examples_dir, "euler_1d_sod", "sod.case"),
        mesh_file=join(examples_dir, "euler_1d_sod", "box.nmsh"),
        user_file=join(examples_dir, "euler_1d_sod", "sod.f90")
    ),
    "euler_2d_cylinder": NekoTestCase(
        case_file=join(examples_dir, "euler_2d_cylinder", "euler_2d_cylinder.case"),
        mesh_file=join(examples_dir, "euler_2d_cylinder", "cyl.nmsh"),
        user_file=join(examples_dir, "euler_2d_cylinder", "euler_2d_cylinder.f90")
    ),
    "euler_2d_forward_facing_step": NekoTestCase(
        case_file=join(examples_dir, "euler_2d_forward_facing_step", "step.case"),
        mesh_file=join(examples_dir, "euler_2d_forward_facing_step", "step.nmsh"),
        user_file=join(examples_dir, "euler_2d_forward_facing_step", "step.f90")
    ),
    "euler_2d_smooth": NekoTestCase(
        case_file=join(examples_dir, "euler_2d_smooth", "euler_2d_smooth.case"),
        mesh_file=join(examples_dir, "euler_2d_smooth", "box.nmsh"),
        user_file=join(examples_dir, "euler_2d_smooth", "euler_2d_smooth.f90")
    ),
    "hemi": NekoTestCase(
        case_file=join(examples_dir, "hemi", "hemi.case"),
        mesh_file=join(examples_dir, "hemi", "hemi.nmsh")
    ),
    "lid": NekoTestCase(
        case_file=join(examples_dir, "lid", "lid.case"),
        mesh_file=join(examples_dir, "lid", "lid.nmsh"),
        user_file=join(examples_dir, "lid", "lid.f90")
    ),
    "rayleigh_benard": NekoTestCase(
        case_file=join(examples_dir, "rayleigh_benard", "rayleigh.case"),
        mesh_file=join(examples_dir, "rayleigh_benard", "box.nmsh"),
        user_file=join(examples_dir, "rayleigh_benard", "rayleigh.f90")
    ),
    "rayleigh_benard_cylinder": NekoTestCase(
        case_file=join(examples_dir, "rayleigh_benard_cylinder", "rayleigh.case"),
        mesh_file=join(examples_dir, "rayleigh_benard_cylinder", "rayleigh.nmsh"),
        user_file=join(examples_dir, "rayleigh_benard_cylinder", "rayleigh.f90")
    ),
    "scalar_mms": NekoTestCase(
        case_file=join(examples_dir, "scalar_mms", "scalar_mms.case"),
        mesh_file=join(examples_dir, "scalar_mms", "box.nmsh"),
        user_file=join(examples_dir, "scalar_mms", "scalar_mms.f90")
    ),
    "tgv": NekoTestCase(
        case_file=join(examples_dir, "tgv", "tgv.case"),
        user_file=join(examples_dir, "tgv", "tgv.f90"),
        mesh_file=join(examples_dir, "tgv", "512.nmsh"),
    ),
    "turb_channel": NekoTestCase(
        case_file=join(examples_dir, "turb_channel", "turb_channel.case"),
        mesh_file=join(examples_dir, "turb_channel", "box.nmsh"),
        user_file=join(examples_dir, "turb_channel", "turb_channel.f90")
    ),
    "turb_pipe": NekoTestCase(
        case_file=join(examples_dir, "turb_pipe", "turb_pipe.case"),
        mesh_file=join(examples_dir, "turb_pipe", "turb_pipe.nmsh"),
        user_file=join(examples_dir, "turb_pipe", "turb_pipe.f90")
    ),
    "programming_custom_types": NekoTestCase(
        user_file=join(examples_dir, "programming", "custom_types.f90")
    ),
    "programming_fields_vectors_math": NekoTestCase(
        user_file=join(examples_dir, "programming", "fields_vectors_math.f90")
    ),
    "programming_output": NekoTestCase(
        user_file=join(examples_dir, "programming", "output.f90")
    ),
    "programming_registries": NekoTestCase(
        user_file=join(examples_dir, "programming", "registries.f90")
    ),
    "programming_startup_json": NekoTestCase(
        user_file=join(examples_dir, "programming", "startup_and_json.f90")
    ),
    "programming_user_file_template": NekoTestCase(
        user_file=join(examples_dir, "programming", "user_file_template.f90")
    ),
}


def manipulate_case(example, case, tmp_path):
    """Change the end_time in the case file to be twice the timestep and set the
    mesh file.

    Parameters
    ----------
    example : str
        The name of the example to modify.
    case : dict
        The dictionary representation of the case file.
    tmp_path : pathlib.Path
        The temporary path where the output files will go.

    """

    case_object = case["case"]
    time_object = case_object["time"]
    timestep = time_object.get("timestep", time_object.get("max_timestep"))
    time_object["max_timestep"] = timestep
    time_object["end_time"] = 2 * timestep
    case_object["mesh_file"] = examples[example].mesh_file
    case_object["output_directory"] = str(tmp_path)


@pytest.mark.parametrize("example", ["hemi", "rayleigh_benard", "cylinder"])
#@pytest.mark.parametrize("example", examples.keys())
def test_example_smoke(example, launcher_script, request, log_file, tmp_path):
    """Run a smoke test for the specified Neko example.

    Parameterized against the examples dictionary keys. Note that not all
    examples will work, because this does not run prepare scripts or copy over
    any extra files like the .bin for the TS_channel example. This can be
    fixed, in principle.

    """
    # Max number of ranks to launch on
    max_nprocs = 2
    nprocs = configure_nprocs(max_nprocs)

    test_name = request.node.name
    case_file = examples[example].case_file
    with open(case_file, "r") as f:
        case = json5.load(f)

    manipulate_case(example, case, tmp_path)

    case_file = join("tests", "test_examples", test_name + ".case")

    with open(case_file, "w") as f:
        json.dump(case, f, indent=4)

    neko = turboneko
    makeneko = get_makeneko()

    if examples[example].user_file:
        result = subprocess.run(
            [makeneko, examples[example].user_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True)
        neko = "./neko"

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"


@pytest.mark.parametrize("example", [key for key in examples.keys() if examples[key].user_file])
def test_example_compile(example, log_file):
    """Compile all examples that have a user file.

    """

    makeneko = get_makeneko()

    with open(log_file, "w") as f:
        result = subprocess.run(
            [makeneko, examples[example].user_file],
            stdout=f,
            stderr=subprocess.STDOUT,
            text=True)
    assert (
        result.returncode == 0
    ), f"makeneko process failed with exit code {result.returncode}"


def test_example_poisson(log_file):
    """The Poisson example is special since it needs to be compiled with make
    and creates its own program. Here, we test compilation.

    """

    result = subprocess.run(
            ["make", "-C", join(examples_dir, "poisson")],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True)

    # Write the output to the log file
    with open(log_file, "w") as f:
        f.write(result.stdout)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"compiling the Poisson example failed with exit code {result.returncode}"
