import os
import subprocess

def get_neko():
    """
    Returns the path to the turboneko executable via the
    environmental variable NEKO_BIN, or just returns "neko" 
    """

    src_dir = os.getenv("NEKO_BIN")
    if src_dir:
        return os.path.join(src_dir, "neko")
    else:
        return "neko"

def get_neko_dir():
    return ".."


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
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)

    return result