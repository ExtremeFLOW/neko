#!/usr/bin/env python3

from __future__ import annotations

import re
import shutil
import sys
import tempfile
from pathlib import Path


CREATE_SUITE_USAGE = """Usage: contrib/add_unit_test/add_unit_test.sh <test_name> <is_parallel>

Creates a new serial or parallel pFUnit test from tests/unit/templates and
wires it into tests/unit/Makefile.am, configure.ac, and tests/unit/.gitignore.

The test name must match: ^[a-z][a-z0-9_]*$
The parallel flag accepts: true/false, yes/no, 1/0
"""

ADD_FILE_USAGE = """Usage: contrib/add_unit_test/add_file_to_unit_test.sh <suite_name> <pf_name>

Creates tests/unit/<suite_name>/test_<pf_name>.pf from the appropriate pFUnit
template and wires it into the existing suite Makefile.in and
tests/unit/Makefile.am.

Both names must match: ^[a-z][a-z0-9_]*$
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


def validate_name(name: str, what: str) -> None:
    """Validate user-provided suite and file names against repo conventions."""
    if not re.fullmatch(r"[a-z][a-z0-9_]*", name):
        die(f"{what} must match ^[a-z][a-z0-9_]*$")


def read_lines(path: Path) -> list[str]:
    """Read a UTF-8 text file and preserve line endings."""
    return path.read_text(encoding="utf-8").splitlines(keepends=True)


def write_text_atomic(path: Path, content: str) -> None:
    """Atomically replace a text file by writing through a sibling temporary file.

    In other words, it writes to a temporary first, so if that fails, the
    original files are never corrupted, and the temporary just gets cleaned up.

    This is probably a bit overkill-safe, but why not.
    """
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", text=True)
    tmp_path = Path(tmp_name)
    try:
        with open(fd, "w", encoding="utf-8") as handle:
            handle.write(content)
        if path.exists():
            tmp_path.chmod(path.stat().st_mode & 0o777)
        else:
            tmp_path.chmod(0o644)
        tmp_path.replace(path)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def write_lines_atomic(path: Path, lines: list[str]) -> None:
    """Atomically replace a text file from an in-memory list of lines."""
    write_text_atomic(path, "".join(lines))


def insert_before_pattern(lines: list[str], pattern: str, new_line: str) -> list[str]:
    """Insert a line before the first matching anchor unless already present."""

    # If the new line is already present, assume the change was made and do not insert
    if new_line in lines:
        return lines

    # Otherwise, find the first line containing the pattern and insert before it.
    for index, line in enumerate(lines):
        if pattern in line:
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

    # If the new line is already present, assume the change was made and do not insert
    if new_line in lines:
        return lines

    # Go through lines, first search for the block_start, and then a line with
    # the patter before the block_end.
    in_block = False
    for index, line in enumerate(lines):
        if block_start in line:
            in_block = True
        if in_block and pattern in line:
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
    """Drop copied build products so a new test starts from a clean template.

    This is run after copying the template directory to create the new test,
    so that the new directory is free of build products in the case when
    the template test was compiled before.
    """
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


def repo_paths() -> tuple[Path, Path]:
    """Return the repository root and the unit-test directory."""
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    unit_dir = repo_root / "tests" / "unit"
    return repo_root, unit_dir


def suite_is_parallel(makefile_lines: list[str]) -> bool:
    """Infer whether an existing suite is MPI-enabled from its Makefile.in."""
    return any(line.strip() == "USEMPI=YES" for line in makefile_lines)


def patch_template_pf(path: Path, test_name: str, is_parallel: bool) -> None:
    """Patch a copied template .pf file to match the requested test name."""
    if is_parallel:
        replace_all(
            path,
            [
                ("module test_parallel", f"module test_{test_name}"),
                ("end module test_parallel", f"end module test_{test_name}"),
                (
                    "type, extends(MPITestCase) :: parallel_template_case",
                    f"type, extends(MPITestCase) :: test_{test_name}_case",
                ),
                ("end type parallel_template_case", f"end type test_{test_name}_case"),
                ("test_parallel_template_passes", "test_template_passes"),
                ("class(parallel_template_case)", f"class(test_{test_name}_case)"),
            ],
        )
        return

    replace_all(
        path,
        [
            ("module test_serial", f"module test_{test_name}"),
            ("end module test_serial", f"end module test_{test_name}"),
            ("test_serial_template_passes", "test_template_passes"),
            ("serial_template_case", f"test_{test_name}_case"),
        ],
    )


def rename_template_files(test_dir: Path, test_name: str, is_parallel: bool) -> tuple[str, str]:
    """Rename template files and patch their contents for the requested suite name."""
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
        patch_template_pf(test_dir / pf_file, test_name, is_parallel=True)
        return test_program, suite_program

    (test_dir / "test_serial.pf").rename(test_dir / pf_file)
    replace_all(
        test_dir / "Makefile.in",
        [
            ("serial_test", test_program),
            ("test_serial.pf", pf_file),
        ],
    )
    patch_template_pf(test_dir / pf_file, test_name, is_parallel=False)
    return test_program, test_program


def add_pf_to_suite_makefile(lines: list[str], pf_file: str) -> list[str]:
    """Append a new .pf entry to a suite-local _TESTS list.

    The function preserves the existing indentation style of the final .pf entry
    so that both tab-aligned and space-aligned lists stay readable.
    """
    if any(pf_file in line for line in lines):
        return lines

    tests_line_index = None
    for index, line in enumerate(lines):
        if "_TESTS :=" in line:
            tests_line_index = index
            break
    if tests_line_index is None:
        raise ValueError("could not find a _TESTS := block in suite Makefile.in")

    last_entry_index = tests_line_index
    for index in range(tests_line_index + 1, len(lines)):
        line = lines[index]
        stripped = line.strip()
        if not stripped:
            break
        if "_OTHER_LIBRARIES" in line or "$(eval $(call make_pfunit_test" in line:
            break
        if ".pf" in line:
            last_entry_index = index
            continue
        if stripped.endswith("\\"):
            last_entry_index = index
            continue
        break

    previous_line = lines[last_entry_index].rstrip("\n")
    indent_match = re.match(r"^(\s*)", previous_line)
    indent = indent_match.group(1) if indent_match else ""
    if last_entry_index == tests_line_index or not indent:
        indent = "    "

    if not previous_line.rstrip().endswith("\\"):
        lines[last_entry_index] = f"{previous_line}\\\n"
    else:
        lines[last_entry_index] = f"{previous_line}\n"

    new_line = f"{indent}{pf_file}\n"
    return lines[: last_entry_index + 1] + [new_line] + lines[last_entry_index + 1 :]


def create_suite(test_name: str, is_parallel: bool) -> None:
    """Create a new unit-test directory and register it in the build system."""
    validate_name(test_name, "test name")
    repo_root, unit_dir = repo_paths()
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

        # Copy template to a temporary directory
        shutil.copytree(template_dir, temp_test_dir)

        # Remove all build artifacts potentially  present
        remove_build_artifacts(temp_test_dir)

        # Rename the template files and patch their contents
        test_program, suite_or_binary = rename_template_files(
            temp_test_dir, test_name, is_parallel
        )


        #
        # Start patching files, everything is done in memory first
        #

        # Patch the Makefile.am under tests/unit
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

        # A little shaky because relies on a comment in the Makefile,
        # maybe change this later.
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

        # Patch configure.ac
        configure_ac = repo_root / "configure.ac"
        configure_original = configure_ac.read_text(encoding="utf-8")
        configure_lines = configure_original.splitlines(keepends=True)
        configure_lines = insert_before_pattern(
            configure_lines,
            "tests/unit/templates/serial/Makefile",
            f"        tests/unit/{test_name}/Makefile\\\n",
        )

        # Patch gitignore
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

        # Write the files that are just lines in memory to their final location
        # via a temporary
        write_lines_atomic(makefile_am, makefile_am_lines)
        write_lines_atomic(configure_ac, configure_lines)
        write_lines_atomic(gitignore, gitignore_lines)
        temp_test_dir.replace(test_dir)
    except Exception:
        # Restore tracked files if any write fails, then remove the temporary
        # directory so the repository stays unchanged on error.
        for path, original in modified_files.items():
            if path.exists():
                write_text_atomic(path, original)
        if temp_test_dir.exists():
            shutil.rmtree(temp_test_dir)
        raise

    print(f"Created tests/unit/{test_name} from templates/{template_name}")


def add_file_to_suite(suite_name: str, pf_name: str) -> None:
    """Create a new .pf file in an existing suite and wire it into the build."""
    validate_name(suite_name, "suite name")
    validate_name(pf_name, "pf name")

    _, unit_dir = repo_paths()
    suite_dir = unit_dir / suite_name
    makefile_in = suite_dir / "Makefile.in"
    pf_file = f"test_{pf_name}.pf"
    pf_path = suite_dir / pf_file
    makefile_am = unit_dir / "Makefile.am"

    if not suite_dir.is_dir():
        die(f"missing suite directory: {suite_dir}")
    if not makefile_in.is_file():
        die(f"missing suite Makefile.in: {makefile_in}")
    if pf_path.exists():
        die(f"target already exists: {pf_path}")

    makefile_lines = read_lines(makefile_in)
    is_parallel = suite_is_parallel(makefile_lines)
    makefile_am_original = makefile_am.read_text(encoding="utf-8")
    makefile_in_original = makefile_in.read_text(encoding="utf-8")
    makefile_am_lines = makefile_am_original.splitlines(keepends=True)
    extra_dist_line = f"\t{suite_name}/{pf_file}\\\n"
    template_dir = unit_dir / "templates" / ("parallel" if is_parallel else "serial")
    template_pf = template_dir / ("test_parallel.pf" if is_parallel else "test_serial.pf")
    if not template_pf.is_file():
        die(f"missing template file: {template_pf}")

    modified_files = {
        makefile_in: makefile_in_original,
        makefile_am: makefile_am_original,
    }

    try:
        # Write the new file first, then update the two tracked build files. If
        # anything fails, roll all three changes back together.
        shutil.copy2(template_pf, pf_path)
        patch_template_pf(pf_path, pf_name, is_parallel)
        updated_makefile_lines = add_pf_to_suite_makefile(makefile_lines, pf_file)
        updated_makefile_am_lines = insert_before_pattern_in_block(
            makefile_am_lines,
            "EXTRA_DIST =",
            "UNIT_TEST_MAKEFILES =",
            "templates/parallel/test_parallel.pf",
            extra_dist_line,
        )

        write_lines_atomic(makefile_in, updated_makefile_lines)
        write_lines_atomic(makefile_am, updated_makefile_am_lines)
    except Exception:
        # Restore tracked files and remove the new .pf on any failure.
        for path, original in modified_files.items():
            if path.exists():
                write_text_atomic(path, original)
        if pf_path.exists():
            pf_path.unlink()
        raise

    print(f"Created tests/unit/{suite_name}/{pf_file}")


def main() -> None:
    """Dispatch between suite creation and file insertion modes."""
    if len(sys.argv) < 2:
        print(CREATE_SUITE_USAGE, end="", file=sys.stderr)
        print(file=sys.stderr)
        print(ADD_FILE_USAGE, end="", file=sys.stderr)
        raise SystemExit(1)

    command = sys.argv[1]

    if command == "create-suite":
        if len(sys.argv) != 4:
            print(CREATE_SUITE_USAGE, end="", file=sys.stderr)
            raise SystemExit(1)
        create_suite(sys.argv[2], parse_bool(sys.argv[3]))
        return

    if command == "add-file":
        if len(sys.argv) != 4:
            print(ADD_FILE_USAGE, end="", file=sys.stderr)
            raise SystemExit(1)
        add_file_to_suite(sys.argv[2], sys.argv[3])
        return

    die(f"unknown command: {command}")


if __name__ == "__main__":
    main()
