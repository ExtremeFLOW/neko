"""Contains useful functions for running neko, makeneko, getting paths, etc.

"""
import os
from os.path import join
import subprocess
from conftest import logger, MAX_NPROCS



def get_neko():
    """
    Returns the path to the turboneko executable via the environmental variable
    NEKO_EXEC, or just returns "neko", assuming it is in the PATH.
    """

    if os.getenv("NEKO_EXEC"):
        return os.getenv("NEKO_EXEC")
    else:
        return "neko"

def get_makeneko():
    if os.getenv("MAKENEKO_EXEC"):
        return os.getenv("MAKENEKO_EXEC")
    else:
        return "makeneko"

def get_neko_dir():
    """Returns the root of the neko dirctory structure relative to the directory
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

    with open(log_file, "w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                                    text=True)

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
