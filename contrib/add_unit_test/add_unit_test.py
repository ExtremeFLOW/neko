#!/usr/bin/env python3

from __future__ import annotations

import re
import shutil
import sys
import tempfile
from pathlib import Path


USAGE = """Usage: contrib/add_unit_test/add_unit_test.sh <test_name> <is_parallel>

Creates a new serial or parallel pFUnit test from tests/unit/templates and
wires it into tests/unit/Makefile.am, configure.ac, and tests/unit/.gitignore.

The test name must match: ^[a-z][a-z0-9_]*$
The parallel flag accepts: true/false, yes/no, 1/0
"""


def die(message: str) -> None:
    """Print an error message and exit with a failure status."""
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(1)


def parse_bool(value: str) -> bool:
    """Parse the user-provided parallel flag into a boolean."""
    lowered = value.lower()
    if lowered in {"true", "yes", "1"}:
        return True
    if lowered in {"false", "no", "0"}:
        return False
    die("parallel flag must be one of: true, false, yes, no, 1, 0")


def read_lines(path: Path) -> list[str]:
    """Read a text file as UTF-8 and preserve line endings."""
    return path.read_text(encoding="utf-8").splitlines(keepends=True)


def write_lines(path: Path, lines: list[str]) -> None:
    """Write lines back to a UTF-8 text file."""
    path.write_text("".join(lines), encoding="utf-8")


def insert_before_pattern(lines: list[str], pattern: str, new_line: str) -> list[str]:
    """Insert a line before the first matching anchor, unless already present there."""
    for index, line in enumerate(lines):
        if pattern in line:
            if index > 0 and lines[index - 1] == new_line:
                return lines
            return lines[:index] + [new_line] + lines[index:]
    raise ValueError(f"could not find anchor: {pattern}")


def insert_before_pattern_in_block(
    lines: list[str],
    block_start: str,
    block_end: str,
    pattern: str,
    new_line: str,
) -> list[str]:
    """Insert a line before an anchor that must appear within a specific block."""
    in_block = False
    for index, line in enumerate(lines):
        if block_start in line:
            in_block = True
        if in_block and pattern in line:
            if index > 0 and lines[index - 1] == new_line:
                return lines
            return lines[:index] + [new_line] + lines[index:]
        if in_block and block_end in line:
            in_block = False
    raise ValueError(f"could not find anchor in block: {pattern}")


def replace_all(path: Path, replacements: list[tuple[str, str]]) -> None:
    """Apply a sequence of plain-text replacements to a UTF-8 file."""
    content = path.read_text(encoding="utf-8")
    for old, new in replacements:
        content = content.replace(old, new)
    path.write_text(content, encoding="utf-8")


def remove_build_artifacts(test_dir: Path) -> None:
    """Drop copied build products so a new test starts from a clean template."""
    patterns = (
        "Makefile",
        "*.o",
        "*.mod",
        "*.a",
        "*.inc",
        "*.F90",
        "*.log",
        "*.trs",
        "serial_test",
        "parallel_suite",
    )
    for pattern in patterns:
        for path in test_dir.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()


def write_text_atomic(path: Path, content: str) -> None:
    """Atomically replace a text file by writing through a sibling temporary file."""
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", text=True)
    tmp_path = Path(tmp_name)
    try:
        with open(fd, "w", encoding="utf-8") as handle:
            handle.write(content)
        tmp_path.replace(path)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def write_lines_atomic(path: Path, lines: list[str]) -> None:
    """Atomically replace a text file from an in-memory list of lines."""
    write_text_atomic(path, "".join(lines))


def rename_template_files(test_dir: Path, test_name: str, is_parallel: bool) -> tuple[str, str]:
    """Rename template files and patch their contents for the requested test name."""
    pf_file = f"test_{test_name}.pf"
    test_program = f"{test_name}_test"

    if is_parallel:
        suite_program = f"{test_name}_suite"
        (test_dir / "test_parallel.pf").rename(test_dir / pf_file)
        (test_dir / "parallel_test").rename(test_dir / test_program)

        replace_all(
            test_dir / "Makefile.in",
            [
                ("parallel_suite", suite_program),
                ("test_parallel.pf", pf_file),
            ],
        )
        replace_all(
            test_dir / test_program,
            [("./templates/parallel/parallel_suite", f"./{test_name}/{suite_program}")],
        )
        replace_all(
            test_dir / pf_file,
            [
                ("module test_parallel", f"module test_{test_name}"),
                ("end module test_parallel", f"end module test_{test_name}"),
                (
                    "type, extends(MPITestCase) :: parallel_template_case",
                    f"type, extends(MPITestCase) :: test_{test_name}_case",
                ),
                (
                    "end type parallel_template_case",
                    f"end type test_{test_name}_case",
                ),
                ("test_parallel_template_passes", "test_template_passes"),
                (
                    "class(parallel_template_case)",
                    f"class(test_{test_name}_case)",
                ),
            ],
        )
        return test_program, suite_program

    (test_dir / "test_serial.pf").rename(test_dir / pf_file)
    replace_all(
        test_dir / "Makefile.in",
        [
            ("serial_test", test_program),
            ("test_serial.pf", pf_file),
        ],
    )
    replace_all(
        test_dir / pf_file,
        [
            ("module test_serial", f"module test_{test_name}"),
            ("end module test_serial", f"end module test_{test_name}"),
            ("test_serial_template_passes", "test_template_passes"),
        ],
    )
    return test_program, test_program


def main() -> None:
    """Create a new unit-test directory and register it in the build system."""
    if len(sys.argv) != 3:
        print(USAGE, end="", file=sys.stderr)
        raise SystemExit(1)

    test_name = sys.argv[1]
    is_parallel = parse_bool(sys.argv[2])

    if not re.fullmatch(r"[a-z][a-z0-9_]*", test_name):
        print(USAGE, end="", file=sys.stderr)
        die("test name must match ^[a-z][a-z0-9_]*$")

    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    unit_dir = repo_root / "tests" / "unit"
    template_name = "parallel" if is_parallel else "serial"
    template_dir = unit_dir / "templates" / template_name
    test_dir = unit_dir / test_name
    temp_test_dir = unit_dir / f".{test_name}.tmp"

    if not template_dir.is_dir():
        die(f"missing template directory: {template_dir}")
    if test_dir.exists():
        die(f"target already exists: {test_dir}")
    if temp_test_dir.exists():
        die(f"temporary target already exists: {temp_test_dir}")

    modified_files: dict[Path, str] = {}

    try:
        # Build the new test in a temporary directory so failures do not leave a
        # partially-created test tree behind.
        shutil.copytree(template_dir, temp_test_dir)
        remove_build_artifacts(temp_test_dir)
        test_program, suite_or_binary = rename_template_files(
            temp_test_dir, test_name, is_parallel
        )

        makefile_am = unit_dir / "Makefile.am"
        makefile_am_original = makefile_am.read_text(encoding="utf-8")
        makefile_am_lines = makefile_am_original.splitlines(keepends=True)
        makefile_am_lines = insert_before_pattern_in_block(
            makefile_am_lines,
            "SUBDIRS =",
            "TESTS =",
            "templates/serial",
            f"\t  {test_name}\\\n",
        )
        makefile_am_lines = insert_before_pattern_in_block(
            makefile_am_lines,
            "TESTS =",
            "# Note we need to list .pf and runner scripts manually",
            "templates/serial/serial_test",
            f"\t{test_name}/{test_program}\\\n",
        )

        extra_dist_anchor = (
            "templates/parallel/test_parallel.pf"
            if is_parallel
            else "templates/serial/test_serial.pf"
        )
        makefile_am_lines = insert_before_pattern_in_block(
            makefile_am_lines,
            "EXTRA_DIST =",
            "UNIT_TEST_MAKEFILES =",
            extra_dist_anchor,
            f"\t{test_name}/test_{test_name}.pf\\\n",
        )
        if is_parallel:
            makefile_am_lines = insert_before_pattern_in_block(
                makefile_am_lines,
                "EXTRA_DIST =",
                "UNIT_TEST_MAKEFILES =",
                extra_dist_anchor,
                f"\t{test_name}/{test_program}\\\n",
            )

        configure_ac = repo_root / "configure.ac"
        configure_original = configure_ac.read_text(encoding="utf-8")
        configure_lines = configure_original.splitlines(keepends=True)
        configure_lines = insert_before_pattern(
            configure_lines,
            "tests/unit/templates/serial/Makefile",
            f"        tests/unit/{test_name}/Makefile\\\n",
        )

        gitignore = unit_dir / ".gitignore"
        gitignore_original = gitignore.read_text(encoding="utf-8")
        gitignore_lines = gitignore_original.splitlines(keepends=True)
        ignore_entry = (
            f"{test_name}/{suite_or_binary}\n"
            if is_parallel
            else f"{test_name}/{test_program}\n"
        )
        if ignore_entry not in gitignore_lines:
            gitignore_lines.append(ignore_entry)

        modified_files = {
            makefile_am: makefile_am_original,
            configure_ac: configure_original,
            gitignore: gitignore_original,
        }

        # Commit tracked-file updates atomically before moving the generated test
        # directory into place.
        write_lines_atomic(makefile_am, makefile_am_lines)
        write_lines_atomic(configure_ac, configure_lines)
        write_lines_atomic(gitignore, gitignore_lines)
        temp_test_dir.replace(test_dir)
    except Exception:
        # Restore the original repository files if any write in the commit phase
        # fails, then remove the temporary test directory.
        for path, original in modified_files.items():
            if path.exists():
                write_text_atomic(path, original)
        if temp_test_dir.exists():
            shutil.rmtree(temp_test_dir)
        raise

    print(f"Created tests/unit/{test_name} from templates/{template_name}")


if __name__ == "__main__":
    main()
