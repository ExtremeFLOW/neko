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

