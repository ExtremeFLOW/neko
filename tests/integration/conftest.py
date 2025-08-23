"""Configuration for pytest integration tests. Contains useful fixtures and
defines command line options.

"""
import os
import pytest
import logging

# The backend used to run Neko.
BACKEND = "cpu"
# Whether the backned is not the CPU.
USES_DEVICE = False
# The max number of ranks to launch on.
MAX_NPROCS = 1e+9
# Real precision
RP = "dp"

logger = logging.getLogger("pytest_configure")
logger.setLevel(logging.DEBUG)  # ensure it captures debug and above

# Check if logger has no handlers and add one
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)  # or .INFO, depending on what you want
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

def pytest_addoption(parser):
    parser.addoption(
        "--launcher-script",
        action="store",
        default="./default_cpu_launcher.sh",  # default value
        help="Path to the launcher script for running neko."
    )
    parser.addoption(
        "--backend",
        action="store",
        default="cpu",
        choices=["cpu", "cuda", "hip", "opencl"],
        help="Backend to use for tests (cpu [default], cuda, hip, opencl)."
    )
    parser.addoption(
        "--max_nprocs",
        action="store",
        default=1e+9,
        help="The maximum number of processes to launch on. Defaults to 1e9 (no limit)."
    )
    parser.addoption(
        "--real_precision",
        action="store",
        default="dp",
        help="Precisions of reals for Neko (dp [default], sp)."
    )

def pytest_configure(config):
    global BACKEND, USES_DEVICE, MAX_NPROCS, RP
    BACKEND = config.getoption("--backend")
    USES_DEVICE = BACKEND != "cpu"
    MAX_NPROCS = int(config.getoption("--max_nprocs"))
    RP = config.getoption("--real_precision")
    print(f"Using backend: {BACKEND}, uses_device: {USES_DEVICE}, max number of ranks: {MAX_NPROCS}, real precision: {RP}")

@pytest.fixture
def launcher_script(request):
    return request.config.getoption("--launcher-script")

@pytest.fixture(scope="session")
def backend(request):
    return request.config.getoption("--backend")

@pytest.fixture(scope="session")
def uses_device(backend):
    return backend != "cpu"

@pytest.fixture
def log_file(request):
    """Fixture to create a log file for each test.

    The log file is named after the test function and stored in a 'logs'
    directory.

    """
    # Use the test function name
    test_name = request.node.name

    # Make sure logs directory exists
    os.makedirs("logs", exist_ok=True)

    # Create the log file path
    log_path = os.path.join("logs", f"{test_name}.log")

    return log_path

