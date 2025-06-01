import os

def pytest_addoption(parser):
    parser.addoption(
        "--launcher-script",
        action="store",
        default="./default_runner.sh",  # default value
        help="Path to the launcher script for running neko."
    )

import pytest

@pytest.fixture
def launcher_script(request):
    return request.config.getoption("--launcher-script")


@pytest.fixture
def log_file(request):
    # Use the test function name
    test_name = request.node.name

    # Make sure logs directory exists
    os.makedirs("logs", exist_ok=True)

    # Create the log file path
    log_path = os.path.join("logs", f"{test_name}.log")

    return log_path
