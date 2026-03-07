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
Neko is simulation software, and a common scenario is that the user will ask
you to help set up a simulation. Alternatively, they may want tips about specific
settings or an explanation about a setup they find in an existing case. The
instructions under this header are specifically for these scenarios. We will
refer to a concrete simulation as "the case".

### General guidelines
Below are the main guidelines you should follow for retrieving information and
making decisions about how to set up a case.

- A simulation case typically lives in its own separate folder, which we will call
  the case folder.
- Input to the simulation consists of at least two files
  - A file with the computational mesh. More instructions on meshing follow.
  - A JSON file, which is either a .json or a .case (most often the latter),
    containing the configuration of the various components of the simulation.
    This file is referred to as the **case file**.
- In many cases, an additional Fortran source file, referred to as the **user
  file** is part of the case setup. This allows advanced customization.
- The folder `examples` contains already set up cases to demonstrate Neko's
  functionality. The `examples/README.md` contains an overview of the cases in
  the different sub-folders, and each of them has a `README.md` of its own. If
  you are tasked with creating a new case, it is a very good idea to base it on
  one of the examples instead of starting from scratch. Think of the examples as
  templates for your work. It is a very good idea to ask the user, whether there
  is a particular example, which would serve as a good starting point for the
  setup.
- When unsure, try to prompt the user for possible details. Most CFD simulations
  will be unique and require a lot of input to setup correctly. "Educated
  guesses" may often not work. If you make assumptions, you should explicitly
  tell the user about them.

### The case file
- Your definitive reference for the contents of the case file is the file
  `doc/pages/user-guide/case-file.md`. You should acknowledge that you got
  access to this file.
- Save the case files as the name of the case folder plus .json, and try to
  pretty-format as good as you can.

### The mesh file
- The general guide for meshing is `doc/pages/user-guide/meshing.md`.
- There are essentially four possible options.
  - The user provides you a mesh file as an .nmsh and you use that.
  - The user has a mesh in .re2 format and you can convert that to .nmsh by
    running `rea2nbin` on that file.
  - The user has a mesh in a different format and you should try to help them
    convert it first to .re2 and then to .nmsh. This is the most difficult
    scenario but you can try to provide guidance based on the meshing.md
    referenced above.
  - There is no mesh, but the geometry is a box. In this case you can generate
    the mesh yourself using the `genmeshbox` utility. The source code for this
    utility is in `contrib/genmeshbox/genmeshbox.f90`. Study it so that you have
    perfect understanding of how it works. Then, follow the user's instructions
    about the mesh: its element sizes (uniform or not), periodicity. Observe
    that for non-uniform element size distributions `genmeshbox` takes the
    location of element edges from files. You may therefore need to generate
    these files based on the user's specifications.
- **Always** run the utility `mesh_checker` on the generated or existing .nmsh
  and make sure it reports no errors.

### The user file
- The general guide for the user file is `doc/pages/user-guide/user-file.md`.
- The file is just a Fortran module.
- Name the user file as the name of the case folder plus .f90.
- If it is possible to set up the case without a user file, you should always do
  that. Only create a user file if it is really necessary for the case setup.
- A **very good** idea is to look into the `examples/programming` directory,
  which contains examples for some basic things you can do in the user file.
  Look at the `README.md` in that subfolder for a description of what is in the
  files, but all the .f90 there are heavily documented.
  Other examples may also contain user files.
- General development guidelines apply:
  - `doc/pages/developer-guide/code-style.md` for the Fortran code style.
  - `doc/pages/developer-guide/important_types.md` provides an overview of
    the most important types for Neko. It is a very good idea to look at their
    implementation in the respective .f90 files under the `src` structure.
- Be generous with comments in the user file you generate, so that one can
  easily follow what your code does.
- At the end always run `makeneko` on the generated user file and make sure it
  compiles.

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
You may be asked to review new code, in the context of a PR or a local git
branch with changes with respect to develop. Your review will be kept to
assessing the criteria outlined below, and not be a general review based on your
opinion about good code. Essentially, you will be an intelligent static code
analyzer, to compensate for the lack of those for Fortran. For the moment, you
will only check Fortran code, and documentation, not CUDA and other accelerator
code. Here is the list of things you need to check.

- For all new files in the PR / branch.
  - Any subroutine that creates a local `allocatable`, must also manually
    `deallocate` it.
  - Any type that has `pointer` or `allocatable` components must have a
    routine, in which the pointers are nullified and the allocatable are
    deallocated. In concrete types, this routine must be called `free`, in types
    with `deferred` procedures, it may be called `free_base` or similar, always
    with a `free_` prefix.
  - If a type has a destructor (`free` or similar) and a constructor (`init` or
    similar), the constructor must begin with calling the destructor.
  - All procedures, including type-bound, must have Doxygen documentation,
    including @param for each dummy argument.
- By default, you should also review documentation. The rules for that are
  found under `doc/AGENTS.md`. If explicitly asked to skip this, follow that
  directive.
