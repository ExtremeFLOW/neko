"""Configuration for pytest integration tests. Contains useful fixtures and 
defines command line options.

"""
import os
import pytest

# The backend used to run Neko.
BACKEND = "cpu"
# Whether the backned is not the CPU.
USES_DEVICE = False

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

def pytest_configure(config):
    global BACKEND, USES_DEVICE
    BACKEND = config.getoption("--backend")
    USES_DEVICE = backend != "cpu"
    print(f"Using backend: {BACKEND}, uses_device: {USES_DEVICE}")

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

