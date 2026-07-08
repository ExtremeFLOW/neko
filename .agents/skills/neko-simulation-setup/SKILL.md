---
name: neko-simulation-setup
description: Help set up, explain, or modify Neko CFD simulation cases, including case files, meshes, and optional user files.
---

# Neko Simulation Setup

Use this skill when helping with a concrete Neko simulation case, advising on
simulation settings, explaining an existing case setup, or creating the files
needed to run a case.

## Core Workflow

- Treat each simulation as a separate case folder.
- The goal is to produce:
  - A computational mesh file, if not already provided by the user.
  - A JSON configuration file, usually `.case` or `.json`, referred to as the
    case file.
  - A Fortran user file, **only** when the setup really needs custom behavior.
- Prefer starting from an existing example instead of creating a case from
  scratch. Read `examples/README.md`, then the relevant example `README.md`.
- Ask the user whether a particular example is a good starting point when this
  is not already clear.
- Prompt for missing CFD-specific details. Avoid relying on educated guesses for
  physical setup, boundary conditions, mesh resolution, or solver choices. When
  making assumptions, state them explicitly.

## Case File

- Use `doc/pages/user-guide/case-file.md` as the definitive reference for case
  file contents, and acknowledge that you got access to this file.
- Name generated case files after the case folder, using `.json` unless the user
  specifically asks for `.case`.
- Pretty-format generated JSON.

## Mesh File

- Read `doc/pages/user-guide/meshing.md` before advising on or creating meshes.
- Handle meshes according to the available source:
  - If the user provides `.nmsh`, use it directly.
  - If the user provides `.re2`, convert it to `.nmsh` with `rea2nbin`.
  - If the user provides another format, help convert it first to `.re2` and
    then to `.nmsh`, using the meshing guide.
  - If there is no mesh and the geometry is a box, generate it with
    `genmeshbox`.
    - Before using `genmeshbox`, study `contrib/genmeshbox/genmeshbox.f90`
      closely enough to understand its inputs and behavior.
    - For non-uniform box meshes, generate edge-location files when needed,
      because `genmeshbox` reads non-uniform element distributions from files.
  - Always run `mesh_checker` on any generated or existing `.nmsh` file and make
  sure it reports no errors.

## User File

- Read `doc/pages/user-guide/user-file.md` before creating or changing a user
  file.
- Create a user file only when the case cannot be expressed cleanly through the
  case file.
- Name the user file after the case folder with a `.f90` extension.
- Look at `examples/programming/README.md` and the documented `.f90` examples
  there before writing custom user-file logic. Other examples may also contain
  relevant user files.
- Follow `doc/pages/developer-guide/code-style.md` for Fortran style.
- Use `doc/pages/developer-guide/important_types.md` and the corresponding
  source files under `src` when working with important Neko types.
- Be generous with explanatory comments in generated user files.
- Run `makeneko` on any generated user file and make sure it compiles.

## Completion Checks

- Confirm the case folder contains the expected case file, mesh file, and any
  necessary user file.
- Confirm generated JSON is formatted and named consistently.
- Confirm `mesh_checker` passed for the `.nmsh` mesh.
- Confirm `makeneko` passed when a user file was created or modified.
- If any check cannot be run locally, clearly tell the user what remains to be
  run and why.
