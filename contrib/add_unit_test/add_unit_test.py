#!/usr/bin/env python3

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


USAGE = """Usage: contrib/add_unit_test/add_unit_test.sh <test_name> <is_parallel>

Creates a new serial or parallel pFUnit test from tests/unit/templates and
wires it into tests/unit/Makefile.am, configure.ac, and tests/unit/.gitignore.

The test name must match: ^[a-z][a-z0-9_]*$
The parallel flag accepts: true/false, yes/no, 1/0
"""


def die(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(1)


def parse_bool(value: str) -> bool:
    lowered = value.lower()
    if lowered in {"true", "yes", "1"}:
        return True
    if lowered in {"false", "no", "0"}:
        return False
    die("parallel flag must be one of: true, false, yes, no, 1, 0")


def read_lines(path: Path) -> list[str]:
    return path.read_text().splitlines(keepends=True)


def write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("".join(lines))


def insert_before_pattern(lines: list[str], pattern: str, new_line: str) -> list[str]:
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
    content = path.read_text()
    for old, new in replacements:
        content = content.replace(old, new)
    path.write_text(content)


def remove_build_artifacts(test_dir: Path) -> None:
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


def rename_template_files(test_dir: Path, test_name: str, is_parallel: bool) -> tuple[str, str]:
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

    if not template_dir.is_dir():
        die(f"missing template directory: {template_dir}")
    if test_dir.exists():
        die(f"target already exists: {test_dir}")

    shutil.copytree(template_dir, test_dir)
    remove_build_artifacts(test_dir)
    test_program, suite_or_binary = rename_template_files(test_dir, test_name, is_parallel)

    makefile_am = unit_dir / "Makefile.am"
    makefile_am_lines = read_lines(makefile_am)
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
    write_lines(makefile_am, makefile_am_lines)

    configure_ac = repo_root / "configure.ac"
    configure_lines = read_lines(configure_ac)
    configure_lines = insert_before_pattern(
        configure_lines,
        "tests/unit/templates/serial/Makefile",
        f"        tests/unit/{test_name}/Makefile\\\n",
    )
    write_lines(configure_ac, configure_lines)

    gitignore = unit_dir / ".gitignore"
    ignore_entry = (
        f"{test_name}/{suite_or_binary}\n"
        if is_parallel
        else f"{test_name}/{test_program}\n"
    )
    gitignore_lines = read_lines(gitignore)
    if ignore_entry not in gitignore_lines:
        gitignore_lines.append(ignore_entry)
        write_lines(gitignore, gitignore_lines)

    print(f"Created tests/unit/{test_name} from templates/{template_name}")


if __name__ == "__main__":
    main()
