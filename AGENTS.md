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

### Unit tests with pFUnit
- When asked to create, extend or fix a unit test, use the project-level skill
  in `.agents/skills/neko-unit-test/SKILL.md`. It covers the templates, the
  scaffolding script in `contrib/add_unit_test`, fixtures, death tests, and
  how to build and run a single suite.
- The human-facing guide is `doc/pages/developer-guide/testing.md`.
- Do not write test code before the scaffolded suite compiles and runs. If you
  cannot run `make` yourself, ask the user to do it.

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
