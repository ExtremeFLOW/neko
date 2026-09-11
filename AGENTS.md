# Neko guide for AI agents

## What is Neko?
Neko is computational fluid dynamics (CFD) software based on the Spectral
Element Method. It primarily targets high-fidelity scale-resolving simulations
of turbulent flows.

## High-level repository structure
- `src`, the majority of the source code.
- `doc`, the Doxygen documentation, with additional pages in markdown format.
- `examples`, curated collection of simulation cases, and user files, showcasing
  capabilities.
- `tests`, unit, integration, and system tests, see "writing tests for Neko"
  below.
- `contrib`, misc scripts and utilities for working with Neko.

## Helping users set up new Neko simulation cases
When the user asks for help setting up, explaining, or modifying a Neko
simulation case, use the project-level skill in
`.agents/skills/neko-simulation-setup/SKILL.md`.

## Writing tests for Neko

### General guidelines

- All the tests are located in the directory `tests`.
- There are three types of tests for Neko, all using different frameworks and
  serving a different purpose.
  - The folder `tests/unit` contains unit tests written using
    [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit). These are run
    by CI for every PR.
  - The folder `tests/integration` contains tests written with pytest. Here,
    python and pytest are used to set up Neko cases, run `neko` and `makeneko` as
    subprocesses and then post-process the results. These are run by CI for
    every PR.
  - The folder `tests/reframe` contains nightly tests that are run on a
    supercomputer via a gitlab pipeline. The tests are written using
    [reframe](https://reframe-hpc.readthedocs.io/en/stable/). These are
    validation tests checking that important cases produce the expected output.
  - *Important*: Since Neko can be run in both double and single precision,
    numerical assertions in tests should adjust their tolerance based on the
    kind of `rp`. The `math` module contains the `NEKO_EPS` parameter that is
    set to the correct machine epsilon for `rp`, thus in Fortran code tolerances
    can be conveniently defined via its value.

### Unit tests with pFUnit.
- A guideline for writing these tests is given in
  `doc/pages/developer-guide/testing.md`. Make sure to read it.
- There are fundamentally two types of tests here: those using MPI to run in
  parallel and those that do not. They differ a bit in file naming and other
  aspects mentioned in the `testing.md` referred to above.
  - A prototype for a test that needs MPI is `tests/unit/field`
  - A prototype for a test that doesn't need MPI is
    `tests/unit/time_based_controller`.
- Make sure you understand how to add the test to the build system, that is
  somewhat tedious. Study `doc/pages/developer-guide/testing.md` carefully to
  that end.
- File/module naming: pFUnit generates a driver that calls
   <pf_basename>_suite(), where <pf_basename> is the .pf filename without
   extension. For a smooth link:
   - Keep the module name inside each .pf file identical to the file’s basename
     (e.g., test_case_file_utils.pf must declare module test_case_file_utils).
  - If you use multiple .pf files in one suite, each file’s module name should
    match its basename and be unique.
- Some tests require a simple `mesh_t` to be constructed programmatically. See
  `subroutine test_field_gen_msh` in `tests/unit/field/test_field_parallel.pf`.
- Note that in JSON routines, the path in the file is separated by periods,
  for example `params.value`. NOT `params/value`.
- If the user asks you to create a new unit test, you should first prepare all
the necessary files and add them the build system. The .pf can just be a dummy
at this point, but make sure to make all the file and module names correct. When
you are done, you should run `make check` in the test folder and make sure it
compiles. If you cannot run this command yourself, you should ask the user to do
it. Only when this succeeds should you start populating the .pf with actual
tests.

### Integration tests with pytest
- These tests are located under `tests/integration`.
- To run and write the tests, pytest is used.
- Each test typically runs one or several neko case configurations, launched as
  subprocesses by pytest, and then uses pytest to check the output correctness.
- Key configuration files are `tests/integration/conftest.py` and
  `tests/integration/testlib.py`. Looking at these, plus existing tests, will
  give you a very good idea of how things work.

## Code review
If asked to review code, follow instructions under
`.github/copilot-instructions.md`.
