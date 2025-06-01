"""Contains useful functions for running neko, makeneko, getting paths, the
backend, etc. 

"""
import os
import subprocess

# The backend used to run Neko.
backend = "cpu"
# Whether the backened is not the CPU.
uses_device = False

def pytest_configure(config):
    global backend, uses_device
    backend = config.getoption("--backend")
    uses_device = backend != "cpu"

def get_neko():
    """
    Returns the path to the turboneko executable via the environmental variable 
    NEKO_BIN, or just returns "neko", assuming it is in the PATH.
    """

    bin_dir = os.getenv("NEKO_BIN")
    if bin_dir:
        return os.path.join(bin_dir, "neko")
    else:
        return "neko"

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

def configure_nprocs(nprocs):
    """
    Validates and adjusts the number of processes (nprocs) based on the system
    and the backend.

    """


    global backend

    real_nprocs = nprocs

    if backend in ("cuda", "hip"):
        # Check available GPUs
        num_gpus = 0

        if backend == "cuda":
            try:
                result = subprocess.run(
                    ["nvidia-smi", "-L"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=5,
                    check=True
                )
                num_gpus = len(result.stdout.decode().strip().split("\n"))
            except Exception as e:
                raise RuntimeError("Could not query CUDA GPUs with nvidia-smi. Error: {}".format(e))

        elif backend == "hip":
            # For HIP, use rocminfo or rocm-smi
            try:
                result = subprocess.run(
                    ["rocminfo"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=5,
                    check=True
                )
                # Count "Agent" entries as GPUs
                num_gpus = sum(1 for line in result.stdout.decode().splitlines() if "Agent" in line)
            except Exception as e:
                raise RuntimeError("Could not query HIP GPUs with rocminfo. Error: {}".format(e))

        if num_gpus < nprocs:
            print(f"[WARNING] Requested {nprocs} processes, but only {num_gpus} GPUs detected for {backend} backend. Adjusting nprocs to {num_gpus}.")
            real_nprocs = num_gpus

    return real_nprocs
