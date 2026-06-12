from pathlib import Path
import subprocess
import sys

import pytest


REPO_ROOT = Path(__file__).resolve().parents[4]
VALIDATOR = REPO_ROOT / "contrib" / "validate_case_schema.py"
EXAMPLES_DIR = REPO_ROOT / "examples"
CASE_FILES = sorted(EXAMPLES_DIR.rglob("*.case")) + sorted(EXAMPLES_DIR.rglob("*.json"))


@pytest.mark.parametrize(
    "case_file",
    CASE_FILES,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_example_case_schema(case_file):
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(case_file)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"schema validation failed for {case_file}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
