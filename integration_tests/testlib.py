import os

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


def run_neko(launcher_script, nprocs, case_file, neko):
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

    Returns
    -------
    result : subprocess.CompletedProcess
        The result of the subprocess run.
    """
    import subprocess

    return subprocess.run(
        [launcher_script, str(nprocs), case_file, neko],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )