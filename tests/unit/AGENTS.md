# Unit Test Guide for AI Agents

This file applies to pFUnit work under `tests/unit`.

Read this together with:
- `doc/pages/developer-guide/testing.md`
- the root `AGENTS.md`

## Default Workflow

Use the helper scripts under `contrib/add_unit_test` by default.

For a brand-new suite, use:

```bash
contrib/add_unit_test/add_unit_test.sh <suite_name> <is_parallel>
```

To add another `.pf` file to an existing suite, use:

```bash
contrib/add_unit_test/add_file_to_unit_test.sh <suite_name> <pf_name>
```

These scripts should be the default path. They create the scaffold from the
generic templates and handle the routine build-system plumbing.

After generating the scaffold:
1. run `make check` in the affected suite directory
2. only then replace the dummy test with real assertions

Use manual build-system edits only when:
- the scripts cannot express the change you need
- you are doing unusual follow-up restructuring after generation

## Examples

Notable working examples in the tree:

- minimal serial: `tests/unit/templates/serial`
- minimal MPI: `tests/unit/templates/parallel`
- serial suite: `tests/unit/time_based_controller`
- MPI suite with runner: `tests/unit/field`
- JSON API usage: `tests/unit/json_utils`
- multi-file serial suite: `tests/unit/time_scheme`
- multi-file MPI suite: `tests/unit/registry`

If a neighboring suite already tests the same subsystem style, copy that before
falling back to a generic template.

## Conventions

Follow the existing naming rules:

- folder: `tests/unit/<name>`
- serial executable: `<name>_test`
- MPI suite target: `<name>_suite`
- MPI runner script: `<name>_test`

Each `.pf` file should:
- be wrapped in a `module`
- use a module name matching the `.pf` basename

Example:
- `test_case_file_utils.pf` should declare `module test_case_file_utils`

Keep test subroutine names at 44 characters or fewer.

## Multi-file Suites

The suite generator creates a single-`.pf` scaffold. If the final suite needs
several `.pf` files, first create the suite with `add_unit_test.sh`, then add
the extra files with `add_file_to_unit_test.sh`.

Use `time_scheme` or `registry` as the model for how to organize a multi-file
suite.

For multi-file suites:
- keep each `.pf` focused on one behavior family
- prefer adding new `.pf` files through `add_file_to_unit_test.sh`
- verify the resulting local `Makefile.in` `_TESTS` list is what you expect

If you add a new runner script for MPI, also add it to `EXTRA_DIST`.

## Writing Good Tests

Test the public contract first. Avoid internal fields unless they are the thing
being tested.

Prefer small, single-purpose tests:
- assert the initial state
- perform one operation
- assert the resulting behavior

Before writing assertions, check whether the test also needs:
- device init/finalize
- a mesh fixture
- a reusable helper

For JSON tests, use `tests/unit/json_utils` as the API example, and remember
that JSON paths use dot notation such as `params.value`.

For type-specialized tests, verify the actual source API in `src/` first. Match
the test data types to the specialized key type, value type, and iterator API.

## Verification

After editing tests, run `make check` in the suite directory.

If you changed a local `Makefile.in` in an already-configured tree, regenerate
the local `Makefile` before assuming the edit took effect.
`tests/unit/Makefile.am` provides:

```bash
make refresh-unit-test-makefiles-local
```

If autotools metadata itself is stale, rerun:

```bash
autoreconf -fi
```

If you only changed comments, a rebuild is usually unnecessary. Say that
explicitly.
