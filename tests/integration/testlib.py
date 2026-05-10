"""Contains useful functions for running neko, makeneko, getting paths, etc.

"""
import os
from os.path import join
import subprocess
from conftest import logger, MAX_NPROCS

def which_command(command):
    """
    Mimics the behavior of the bash command "which [cmd]".
    Returns the absolute path of the executable, based on the PATH env.
    variable. If not found, returns None
    """

    path_env = os.getenv("PATH")
    if not path_env:
        return None

    for path in path_env.split(":"):
        test_path = os.path.join(path, command)
        if os.path.isfile(test_path):
            return test_path

    return None

def get_genmeshbox():
    """
    Returns the path to the genmeshbox executable, via the environmental 
    variable GENMESHBOX_EXEC, or returns the path to the first genmeshbox
    executable in the PATH environment variable.
    """
    if os.getenv("GENMESHBOX_EXEC"):
        return os.getenv("GENMESHBOX_EXEC")
    else:
        cmd = which_command("genmeshbox")
        if cmd is None:
            raise FileNotFoundError("The genmeshbox executable could not be found")
        return cmd

def get_neko():
    """
    Returns the path to the turboneko executable via the environmental variable
    NEKO_EXEC, or just returns the path to the first "neko" executable in the 
    PATH environmental variable.
    """
    if os.getenv("NEKO_EXEC"):
        return os.getenv("NEKO_EXEC")
    else:
        cmd = which_command("neko")
        if cmd is None:
            raise FileNotFoundError("The neko executable could not be found")
        return cmd

def get_makeneko():
    """
    Returns the path to the makeneko executable via the environmental variable
    NEKO_EXEC, or just returns the path to the first "makeneko" executable in the 
    PATH environmental variable.
    """
    if os.getenv("MAKENEKO_EXEC"):
        return os.getenv("MAKENEKO_EXEC")
    else:
        cmd = which_command("makeneko")
        if cmd is None:
            raise FileNotFoundError("The makeneko executable could not be found")
        return cmd

def get_neko_dir():
    """Returns the root of the neko directory structure relative to the directory
    with the integration tests

    """
    return "../.."


def run_neko(launcher_script, nprocs, case_file, neko, log_file):
    """
    Runs the neko executable with the specified parameters.

    Parameters
    ----------
    launcher_script : str
        The path to the launcher script.
    nprocs : int
        The number of processes to launch.
    case_file : str
        The path to the case file.
    neko : str
        The path to the neko executable.
    log_file : str,
        Output is written to this file instead of being captured.

    Returns
    -------
    result : subprocess.CompletedProcess
        The result of the subprocess run.
    """
    cmd = [launcher_script, str(nprocs), case_file, neko]

    try:
        with open(log_file, "w") as f:
            result = subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True
            )
    except Exception as e:
        with open(log_file, "a") as f:
            f.write(f"\n[pytest wrapper] Failed to start process:\n{repr(e)}\n")
        raise
    return result

def parse_log(log_file, output_file):
    """
    Parses the neko log file and saves the output to a specified file.

    Parameters
    ----------
    log_file : str
        The path to the neko log file.
    output_file : str
        The path to the output file where parsed data will be saved.
    """

    neko_dir = get_neko_dir()

    command = [
        "python",
         join(neko_dir, "contrib/neko_log_parser/neko_log_parser.py"),
        log_file,
        "-o",
        output_file
    ]

    result = subprocess.run(command, capture_output=True, text=True)

def configure_nprocs(nprocs):
    """
    Validates and adjusts the number of processes (nprocs) based on the CLI
    option and the test.

    Returns min(nprocs, MAX_NPROCS)
    """

    global MAX_NPROCS

    real_nprocs = nprocs
    if nprocs > MAX_NPROCS:
        logger.warning(f"Requested {nprocs} processes is larger than {MAX_NPROCS}. Setting nprocs to {MAX_NPROCS}.")
        real_nprocs = MAX_NPROCS

    return real_nprocs
