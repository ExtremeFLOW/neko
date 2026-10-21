---
name: neko-unit-test
description: Create, extend, build and run pFUnit unit tests for Neko, including MPI suites, mesh fixtures and death tests that assert on neko_error.
---

# Neko Unit Tests with pFUnit

Use this skill whenever the task is to add a unit-test suite under
`tests/unit`, add a `.pf` file to an existing suite, wire a test into the
build system, or write test code that exercises Neko types with pFUnit.
Integration tests (`tests/integration`, pytest) and ReFrame tests
(`tests/reframe`) are out of scope.

## Ground truth

Read these before changing anything. They are the reference for the
current best-practice setup; if this skill and the files disagree, the
files win and this skill should be updated.

- `tests/unit/templates/serial/` and `tests/unit/templates/parallel/`:
  the canonical `Makefile.in`, `.pf` and (for MPI) runner script. Every new
  suite is a copy of one of these.
- `contrib/add_unit_test/`: `add_unit_test.sh` scaffolds a suite from a
  template, `add_pf_to_unit_test.sh` adds a `.pf` file to an existing suite.
  Both drive `add_unit_test.py`, which needs Python 3.7 or newer.
- `tests/unit/fixtures/`: shared helper modules compiled into suites that
  need them. `mesh_fixture.f90` builds small `mesh_t` objects
  programmatically; `error_redirection.f90` turns `neko_error` and
  `neko_warning` into pFUnit exceptions.
- `tests/unit/errors/`: the reference death-test suite.
- `tests/unit/bc_projector/`: the reference for a suite that compiles a
  fixture, and `tests/unit/wall_sampling/` for a suite with several `.pf`
  files and per-test docstrings.
- `src/comm/pfunit_comm_utils.F90`: `comm_init_test` and `comm_free_test`,
  the only supported way to set up Neko's communicator state in MPI suites.
- `doc/pages/developer-guide/testing.md`: the human-facing guide.

## Naming rules

- Suite name: `^[a-z][a-z0-9_]*$`, one directory per class or closely
  related family of classes (see `stack`). The directory does not mirror
  `src`.
- Serial suite: directory `tests/unit/<name>`, binary `<name>_test`.
- MPI suite: binary `<name>_suite` plus a runner script `<name>_test` that
  launches it with `mpirun`/`mpiexec`. Automake runs the runner, so the
  `TESTS` entry is `<name>/<name>_test` for both kinds.
- Test files are `test_<something>.pf` and the module inside must have the
  same name as the file basename, unique within the suite. pFUnit generates a
  driver that calls `<basename>_suite()`, so a mismatch breaks linking.
- Test subroutine names start with `test_` and describe the behaviour, for
  example `test_base_initialises_resolved_state`. Keep them at most 44
  characters; longer names make the preprocessor fail with the unhelpful
  `Syntax error in argument list`.
- Give every test a `!>` docstring stating what it verifies.

## Workflow for a new suite

Work in this order. Do not write real test code before the dummy suite
compiles and runs.

1. Decide whether the code under test needs MPI. Anything touching
   `mesh_t`, `dofmap_t`, `gs_t`, fields, or `NEKO_COMM` is an MPI suite.
   Pure containers, JSON helpers and time schemes are serial.
2. Scaffold from the top of the repository:
   ```sh
   contrib/add_unit_test/add_unit_test.sh <name> <true|false>   # true = MPI
   contrib/add_unit_test/add_pf_to_unit_test.sh <name> <pf_name> # optional extra files
   ```
   The script copies the template, renames files and identifiers, and edits
   four tracked files: `tests/unit/Makefile.am` (`SUBDIRS`, `TESTS`,
   `EXTRA_DIST`), `configure.ac` (`AC_CONFIG_FILES` under `# Config
   tests/unit`), and `tests/unit/.gitignore`. Review the diff with
   `git status` and `git diff`, but do not redo its work by hand.
   If the default `python3` is older than 3.7 the script fails with
   `future feature annotations is not defined`; use a newer interpreter.
3. Regenerate the build. Neko uses `AM_MAINTAINER_MODE`, so a plain `make`
   does not notice the changed `configure.ac` and `Makefile.am`. Run
   ```sh
   ./regen.sh && ./config.status --recheck && ./config.status
   ```
   which re-runs `configure` with the recorded flags. `configure` needs the
   same environment as the original run (compilers, `nvcc`, runtime library
   paths), so if it fails on a dependency check, ask the user to run it.
   Afterwards `tests/unit/<name>/Makefile` must exist.
4. Compile only the new suite. This requires `libneko` to be built already
   (`make -j` at the top level).
   ```sh
   make -C tests/unit/<name> check
   ```
   Note that `check` inside a suite directory only builds the binary; it
   does not run it.
5. Run only the new suite through the automake test driver:
   ```sh
   make -C tests/unit check-TESTS TESTS=<name>/<name>_test
   ```
   The result is printed as `PASS`/`FAIL` and the output goes to
   `tests/unit/<name>/<name>_test.log`. Use `check-TESTS`, not `check`:
   `make check` in `tests/unit` recurses into every suite first.
   Alternatives: a serial binary can be run directly,
   `./tests/unit/<name>/<name>_test -v -f <pattern>` (the driver accepts
   `-f` to filter tests by name and `-v` for verbose output); an MPI suite
   is run from `tests/unit` with `sh <name>/<name>_test`, because the runner
   uses paths relative to that directory. A `cannot open shared object
   file` error means the runtime library path of the configured
   dependencies (json-fortran, HDF5, MPI) is not set in the current shell.
   If you cannot build or run, stop and ask the user to do it before
   writing tests.
6. Replace the template test in the `.pf` with real tests, rebuild, rerun.
   Finish with a full `make check` from the top level when feasible.

## Anatomy of the Makefile.in

The template `Makefile.in` already contains everything a plain suite needs:
the `USEMPI=YES` line (MPI only), the `NEKO_LIB` dependency so test
objects are rebuilt when `libneko` changes, a rule that regenerates the
`Makefile` from `Makefile.in` through `config.status`, the pFUnit
`make_pfunit_test` call, and a `clean` target. Only touch these parts:

- `<binary>_TESTS`: the list of `.pf` files, one per line with `\`
  continuations. `add_pf_to_unit_test.sh` maintains it.
- `<binary>_OTHER_SOURCES`: fixtures compiled into the suite, see below.
- `clean`: add any files the tests write (`*.nmsh`, `*.chkp`, `*.h5`, ...)
  and add matching patterns to `tests/unit/.gitignore`.

Do not copy `Makefile.in` from an arbitrary older suite; several still lack
the `NEKO_LIB` dependency and the `Makefile` regeneration rule.

## Anatomy of the .pf

Serial suites extend `TestCase`; MPI suites extend `MPITestCase`. Keep the
template `setUp`/`tearDown` and extend them, do not replace them:

- Serial: `device_init`/`device_finalize` guarded by
  `NEKO_BCKND_DEVICE .eq. 1`.
- MPI: `comm_init_test(this%getMpiCommunicator())`, then device init, then
  `neko_mpi_types_init()`; `tearDown` undoes them in reverse and ends with
  `comm_free_test()`. Never duplicate the communicator by hand or assign
  `pe_rank`/`pe_size` directly. Suites without a test-case type (`math`,
  `mesh`, `crystal_router`) take `class(MpiTestMethod)` arguments and call
  `comm_init_test` at the start and `comm_free_test` at the end of every
  test instead.
- MPI tests carry `@test(npes=[...])`. The `-np` in the runner script must be
  at least the largest `npes` in the suite, otherwise pFUnit fails the test
  with `Insufficient processes`. Most suites use 1 or 2 ranks; a suite that
  needs more than a typical machine has cores must let MPI oversubscribe,
  see `tests/unit/crystal_router/crystal_router_test`. A runner may also pin
  environment variables when reproducibility needs it, see
  `tests/unit/restart_consistency/restart_consistency_test`.
- Set up expensive objects (mesh, space, dofmap, coefficients) once as
  components of the test-case type in `setUp` and free them in `tearDown`,
  as in `tests/unit/bc/test_bc.pf`.

Numerics and precision:

- Neko builds and tests in both single and double precision. Use
  `real(kind=rp)` and `_rp` literals, and derive tolerances from `NEKO_EPS`
  in the `math` module instead of hard-coded numbers, for example
  `@assertEqual(expected, actual, tolerance=10*NEKO_EPS)` or
  `@assertRelativelyEqual`. Time is `dp` regardless of `rp`.
- On device builds the data lives on the device. Guard host/device copies
  with `NEKO_BCKND_DEVICE .eq. 1` and copy with
  `copy_from(HOST_TO_DEVICE, sync = .true.)` before asserting on host
  arrays, as in `tests/unit/wall_sampling/test_base.pf`.
- JSON paths passed to `json_utils` routines use periods, `params.value`,
  not slashes.

## Meshes

Do not read mesh files in unit tests. Use `tests/unit/fixtures/mesh_fixture.f90`
(`single_unit_hex_mesh`, `single_reference_element_mesh`,
`single_skewed_element_mesh`, `two_adjacent_unit_hex_mesh`,
`three_hex_right_triangle_tip_mesh`) and extend it when a new shape is
needed. Wire it into the suite as `tests/unit/bc_projector/Makefile.in`
does:

```make
../fixtures/mesh_fixture.o : ../fixtures/mesh_fixture.f90 $(NEKO_LIB)
	$(FC) -c $(FFLAGS) -o $@ $<

<binary>_OTHER_SOURCES = ../fixtures/mesh_fixture.f90
test_<x>.o: ../fixtures/mesh_fixture.o      # one line per .pf that uses it
```

and add `../fixtures/mesh_fixture.o` to the `clean` target. Fixture sources
must also be listed in `EXTRA_DIST` in `tests/unit/Makefile.am`.

## Death tests: asserting that neko_error is raised

`neko_error`, `neko_type_error` and `neko_type_registration_error` end by
calling the procedure pointer `throw_error` in `utils`, and `neko_warning`
calls `throw_warning`. Outside tests the first stops the program and the
second only prints. The fixture `tests/unit/fixtures/error_redirection.f90`
can point either at a routine that raises a pFUnit exception instead, so a
test can check that an error or warning was emitted. Errors and warnings are
switched independently: `redirect_errors`/`restore_errors` and
`redirect_warnings`/`restore_warnings`.

Enable error redirection for a whole suite in its `Makefile.in`, exactly as
`tests/unit/errors/Makefile.in` does:

```make
../fixtures/error_redirection.o : ../fixtures/error_redirection.f90 $(NEKO_LIB)
	$(FC) -c $(FFLAGS) -o $@ $<

<binary>_OTHER_SOURCES = ../fixtures/error_redirection.f90
<binary>_EXTRA_USE := error_redirection
<binary>_EXTRA_INITIALIZE := redirect_errors
$(eval $(call make_pfunit_test,<binary>))
<binary>_driver.o: ../fixtures/error_redirection.o
test_<x>.o: ../fixtures/error_redirection.o
```

plus `../fixtures/error_redirection.o` in `clean`. pFUnit calls
`redirect_errors` once when the driver starts, before any test.
Alternatively, drop the two `EXTRA_` lines and call `redirect_errors` in
`setUp` and `restore_errors` in `tearDown`, which lets a death test live in
an ordinary suite.

In the `.pf`, call the routine that must fail and consume the exception
immediately:

```fortran
call neko_error("boom")
@assertExceptionRaised()
```

Warnings are common on valid code paths, so never redirect them for a whole
suite. A test that expects a warning switches the redirection on and off
itself:

```fortran
call redirect_warnings()
call routine_that_warns()
@assertExceptionRaised()
call restore_warnings()
```

Rules that follow from how the redirection works:

- While errors are redirected, any error a test does not consume with
  `@assertExceptionRaised()` fails that test.
- After the redirected `throw`, control returns to the library routine,
  which continues past the `neko_error` call. Only test error paths where
  the routine can safely run to completion afterwards, and do not use its
  results.
- The error text is not propagated to the exception (the `utils` routines
  pass an empty message), so use the bare `@assertExceptionRaised()` and do
  not assert on the message.
- The only existing example, `tests/unit/errors`, is serial. The same
  mechanism applies to MPI suites but has not been exercised there.

## Optional dependencies

A suite that needs HDF5 (or another optional dependency) cannot be
scaffolded by the script alone. Create it with the script, then move its
`SUBDIRS`, `TESTS` and `EXTRA_DIST` entries in `tests/unit/Makefile.am` into
the `if ENABLE_HDF5` block at the end of the file, following `vtkhdf` and
`restart_consistency_hdf5`. The `configure.ac` entry stays unconditional.

## Checklist before handing over

- Module name equals the `.pf` basename in every file of the suite.
- `tests/unit/Makefile.am` lists the directory in `SUBDIRS`, the runner or
  binary in `TESTS`, and every `.pf`, runner script and fixture in
  `EXTRA_DIST`.
- `configure.ac` lists `tests/unit/<name>/Makefile`.
- `tests/unit/.gitignore` lists the binary and any files the tests write;
  `clean` removes the same files.
- `make -C tests/unit/<name> check` builds and
  `make -C tests/unit check-TESTS TESTS=<name>/<name>_test` passes, ideally
  in both `--enable-real=sp` and `dp` builds if a change is precision
  sensitive.
- Tolerances are expressed through `NEKO_EPS`, not literal constants.
- If anything could not be run locally, say exactly which command remains
  to be run and why.
