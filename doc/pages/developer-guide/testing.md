# Testing {#testing}

\tableofcontents

Neko has three kinds of tests, all located under `tests`:

- Unit tests written with [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit)
  in `tests/unit`. CI runs them for every pull request.
- Integration tests written with pytest in `tests/integration`. They set up
  small cases, run `neko` and `makeneko` as subprocesses and check the output.
  CI runs them for every pull request.
- Validation tests written with [ReFrame](https://reframe-hpc.readthedocs.io)
  in `tests/reframe`, run nightly on a supercomputer.

This page covers the unit tests. The AI-agent skill in
`.agents/skills/neko-unit-test/SKILL.md` describes the same procedure for
agents, so you can let an agent scaffold a suite and then fill in the tests
yourself.

## Installing pFUnit {#testing-pfunit-install}

Neko's CI builds against pFUnit v4.15.0; any recent 4.x release should work.
Replace the installation path in the commands below.
```bash
git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git -b v4.15.0
cd pFUnit && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/pfunit_install_path .. && make -j$(nproc) && make install
```
This installs into a directory named after the major and minor version, in
this case `/pfunit_install_path/PFUNIT-4.15`. pFUnit needs a Python 3
interpreter, since the `.pf` files are preprocessed with a Python script.

## Configuring Neko {#testing-configure}

Pass the pFUnit directory to `configure`, for example
```bash
./configure FC=${FC} FCFLAGS="-O2 -pedantic -std=f2008" --with-pfunit=/pfunit_install_path/PFUNIT-4.15
```
The configuration summary printed at the end should say `pFUnit ... yes`.

## Running the tests {#testing-run}

From the top-level directory, run
```bash
make check
```
This builds Neko and all the unit tests, runs them and prints a summary. The
output of each suite is stored in `tests/unit/<suite>/<suite>_test.log`, and
`tests/unit/test-suite.log` collects the failures.

### Running a single suite {#testing-run-single}

Each suite lives in its own directory under `tests/unit` and can be built and
run on its own, provided the Neko library is already built:
```bash
make -C tests/unit/<suite> check                           # only builds the suite
make -C tests/unit check-TESTS TESTS=<suite>/<suite>_test  # only runs the suite
```
Use `check-TESTS` rather than `check` for the second step, since `make check`
in `tests/unit` first recurses into every suite directory.

A serial suite is a plain executable and can also be run directly. The pFUnit
driver accepts `-f` to run only the tests whose names match a pattern and `-v`
for verbose output:
```bash
./tests/unit/<suite>/<suite>_test -v -f <pattern>
```
An MPI suite is launched through its runner script, which uses paths relative
to `tests/unit`:
```bash
cd tests/unit && sh <suite>/<suite>_test
```
If a test executable fails to start with `cannot open shared object file`,
the runtime library path for the dependencies given to `configure` (for
example json-fortran, HDF5 or MPI) is not set in your shell.

## Adding a new suite {#testing-add-suite}

Unlike `src`, the `tests/unit` directory is not organised by source folder.
Each suite directory loosely corresponds to a single class or to a small
family of related classes, see for example `stack`.

Suite names must match `^[a-z][a-z0-9_]*$`. Two flavours exist:

- A serial suite, for code that does not need MPI. The suite is compiled into
  the executable `<suite>_test`.
- An MPI suite, for anything that touches `mesh_t`, `dofmap_t`, `gs_t`,
  fields or `NEKO_COMM`. The suite is compiled into `<suite>_suite` and run
  through the shell script `<suite>_test`, which launches it with `mpirun` or
  `mpiexec`.

Both flavours have a template under `tests/unit/templates`, and the script
`contrib/add_unit_test/add_unit_test.sh` creates a suite from the template and
registers it in the build system. It needs Python 3.7 or newer.
```bash
contrib/add_unit_test/add_unit_test.sh <suite> false   # serial suite
contrib/add_unit_test/add_unit_test.sh <suite> true    # MPI suite
```
The script creates `tests/unit/<suite>` with a `Makefile.in`, a
`test_<suite>.pf` containing one passing test and, for MPI suites, the runner
script. It then adds the suite to `tests/unit/Makefile.am`, `configure.ac` and
`tests/unit/.gitignore`. Review the result with `git status` and `git diff`.

Since `configure.ac` and `Makefile.am` changed, the build system must be
regenerated. 
Run
```bash
./regen.sh && ./config.status --recheck && ./config.status
```
which re-runs `configure` with the flags recorded in `config.status`, or run
`./regen.sh` followed by `./configure` with your usual flags.

Now build and run the new suite as described in
[Running a single suite](@ref testing-run-single). Only once the template test
passes should you replace it with real tests.

### Adding a test file to an existing suite {#testing-add-file}

A suite can contain several `.pf` files, each holding one module. To add one,
run
```bash
contrib/add_unit_test/add_pf_to_unit_test.sh <suite> <name>
```
which creates `tests/unit/<suite>/test_<name>.pf` from the matching template
and adds it to the suite's `Makefile.in` and to `EXTRA_DIST` in
`tests/unit/Makefile.am`.

## Writing tests {#testing-writing}

It is up to you to learn pFUnit from its own documentation. The following
conventions are specific to Neko.

- The module name inside a `.pf` file must equal the file's base name:
  `test_a.pf` contains `module test_a`. pFUnit generates a driver that calls
  `test_a_suite()`, so a mismatch breaks the link.
- Test subroutine names start with `test_` and describe the behaviour under
  test. Keep them at most 44 characters long; longer names make the
  preprocessor fail with the obscure message
  `Error: Syntax error in argument list at (1)`.
- Give every test a `!>` docstring stating what it verifies.
- Keep the `setUp` and `tearDown` routines from the template and extend them
  rather than replacing them.
- Expensive objects that several tests share, such as a mesh, function space,
  dofmap and coefficients, belong as components of the test-case type and are
  built in `setUp` and freed in `tearDown`, see `tests/unit/bc/test_bc.pf`.
- Paths passed to the `json_utils` routines use periods as separators, for
  example `case.fluid.scheme`.

### Serial and MPI suites {#testing-mpi}

A serial test case extends pFUnit's `TestCase`. Its `setUp` calls
`device_init` and its `tearDown` calls `device_finalize`, both guarded by
`NEKO_BCKND_DEVICE .eq. 1`, so that the suite also works in accelerator
builds.

An MPI test case extends `MPITestCase`. Neko's communicator state is set up
with `comm_init_test(this%getMpiCommunicator())` from the module
`pfunit_comm_utils` in `setUp`, followed by the device initialisation and
`neko_mpi_types_init()`. The `tearDown` undoes this in reverse order and ends
with `comm_free_test()`. Never duplicate the communicator by hand or assign
`pe_rank` and `pe_size` directly. Suites that do not define a test-case type
use `class(MpiTestMethod)` arguments and call `comm_init_test` at the
beginning and `comm_free_test` at the end of each test, see
`tests/unit/math/test_math_parallel.pf`.

Each MPI test declares how many ranks it runs on with `@test(npes=[...])`.
The `-np` in the runner script must be at least the largest value used in the
suite, otherwise pFUnit reports `Insufficient processes to run this test`.
Most suites use one or two ranks. A suite that needs more ranks than a typical
machine has cores must let MPI oversubscribe, see
`tests/unit/crystal_router/crystal_router_test`. Runner scripts may also pin
environment variables when a test needs reproducible behaviour, see
`tests/unit/restart_consistency/restart_consistency_test`.

### Precision {#testing-precision}

Neko builds and tests in both single and double precision. Declare reals as
`real(kind=rp)` with `_rp` literals and derive tolerances from `NEKO_EPS` in
the `math` module, which is the machine epsilon of `rp`, rather than from
hard-coded numbers:
```fortran
@assertEqual(expected, actual, tolerance = 100.0_rp * NEKO_EPS)
```
Time is always double precision (`dp`), independently of `rp`.

In accelerator builds the data of Neko's containers lives on the device. Copy
it to the host before asserting on host arrays, guarded by
`NEKO_BCKND_DEVICE .eq. 1`, see `tests/unit/wall_sampling/test_base.pf`.

### Mesh fixtures {#testing-fixtures}

Unit tests should not read mesh files. The module in
`tests/unit/fixtures/mesh_fixture.f90` builds small meshes programmatically:
a single unit hexahedron, the reference element, a skewed element, two
adjacent hexahedra and a three-element mesh with a right-angled tip. Extend it
when a new shape is needed. A suite that uses it compiles the fixture into
the suite through `_OTHER_SOURCES` and explicit object dependencies in its
`Makefile.in`; copy the relevant lines from
`tests/unit/bc_projector/Makefile.in`. Fixture sources must be listed in
`EXTRA_DIST` in `tests/unit/Makefile.am`.

### Testing error handling {#testing-errors}

`neko_error`, `neko_type_error` and `neko_type_registration_error` end by
calling the procedure pointer `throw_error` in the `utils` module, which in a
simulation stops the program, and `neko_warning` calls `throw_warning`, which
only prints. The fixture `tests/unit/fixtures/error_redirection.f90` can point
either at a routine that raises a pFUnit exception instead, so that a test can
check that an error is emitted:
```fortran
call neko_error("boom")
@assertExceptionRaised()
```
Error redirection is enabled for a whole suite in its `Makefile.in` through
pFUnit's `_EXTRA_USE` and `_EXTRA_INITIALIZE` hooks, which call the fixture's
`redirect_errors` before the first test; copy the lines from
`tests/unit/errors/Makefile.in`. Alternatively, call `redirect_errors` in
`setUp` and `restore_errors` in `tearDown`, which lets such tests live in an
ordinary suite.

Warnings are common on valid code paths, so they are never redirected for a
whole suite. A test that expects a warning switches the redirection on and off
itself with `redirect_warnings` and `restore_warnings` around the call.

Keep the following in mind:

- While errors are redirected, every error must be consumed with
  `@assertExceptionRaised()`, otherwise the test fails.
- After the exception is raised, control returns to the library routine, which
  continues past the `neko_error` call. Only test error paths where the
  routine can safely run to completion, and do not use its results.
- The error message is not propagated to the exception, so assert only that
  an exception was raised, not on its text.

### Suites that depend on optional libraries {#testing-optional}

A suite that needs an optional dependency such as HDF5 is created with the
script like any other, after which its `SUBDIRS`, `TESTS` and `EXTRA_DIST`
entries in `tests/unit/Makefile.am` are moved into the corresponding
conditional block at the end of that file, see `if ENABLE_HDF5`. The entry in
`configure.ac` stays unconditional.

## What the script does {#testing-manual}

For reference, or if you need to wire a suite by hand, these are the places a
suite is registered in:

1. `tests/unit/<suite>/Makefile.in`, copied from the template. The
   `check` target, the `_TESTS` and `_OTHER_LIBRARIES` variables, the
   `make_pfunit_test` call and the `clean` target all use the executable name,
   `<suite>_test` for a serial suite and `<suite>_suite` for an MPI suite.
   `configure` turns `Makefile.in` into `Makefile`.
2. For an MPI suite, the runner script `tests/unit/<suite>/<suite>_test`,
   which launches `./<suite>/<suite>_suite` with `mpirun` or `mpiexec`.
3. `tests/unit/Makefile.am`: the directory in `SUBDIRS`, `<suite>/<suite>_test`
   in `TESTS`, and every `.pf` file, runner script and fixture source in
   `EXTRA_DIST`.
4. `configure.ac`: `tests/unit/<suite>/Makefile` in the `AC_CONFIG_FILES` list
   below the comment `# Config tests/unit`.
5. `tests/unit/.gitignore`: the executable, and any files the tests write. The
   suite's `clean` target should remove the same files.
