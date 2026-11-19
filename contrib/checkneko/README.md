## `checkneko`

`checkneko` checks a Neko case file without running Neko, so that mistakes
are found before a simulation starts rather than after it has read the mesh.
It needs neither MPI nor the Neko library.

## Usage

```text
checkneko [-h] [-d dir] [-s dir] file ...
```

- `-d dir` sets the directory Neko will run in, which relative file names in
  the case are relative to. The default is the current directory, as for
  Neko.
- `-s dir` reads the schemas from a directory, such as `doc/schemas`, instead
  of the installed ones, see [Schemas](#schemas).
- `-h` prints a summary of the options.

The options come before the files, and `--` ends them. The exit status is 1
if a case file has errors, and 0 otherwise, also when there are warnings.

## What is checked

1. The case file is parsed with json-fortran, with Neko's settings, so it
   accepts exactly what Neko accepts, including comments and trailing commas.
2. It is validated against the JSON schema of the case file in
   `doc/schemas`, which is installed with the tool. Neko's distinction between
   integers and reals is enforced: write `1.0` where Neko reads a real and
   `1` where it reads an integer. Misspelled keys and values get a
   suggestion.
3. Values that name a constant must name one in `case.constants`, of the
   kind Neko reads there.
4. The files the case reads must exist: the mesh, restart files, field files
   of initial and boundary conditions, probe points and the like.
5. The labeled zones of an `.nmsh` mesh are read, without the elements, and
   the `zone_indices` of the boundary conditions of the fluid and the scalars
   are checked against them the way Neko does: each zone must be in the mesh,
   have one boundary condition and not two. The compressible solver only
   requires the zones to be in the mesh. Two scalars must not have the same
   name.
6. Point zones that values refer to must be defined in `case.point_zones`.
7. Keys given twice, of which Neko uses the first, and unknown keys close to
   a known one give warnings.

A table lists the labeled zones with their number of facets and the types of
the boundary conditions of the fluid and of each scalar on them, with the
`name` of an entry in parentheses, for example

```text
 Zone  Facets  Fluid               Scalar temperature
    1      72  velocity_value      dirichlet
    2      72  outflow+dong        neumann
    7      96  no_slip (cylinder)  dirichlet
```

A field without a boundary condition on a zone shows `-`. A zone that a
boundary condition refers to but the mesh does not have is listed with `-`
facets, and the facets are `?` if the mesh was not read.

Problems are printed after the table as `file:line: error: path: message`,
with the path in the case file such as
`case.fluid.boundary_conditions[2].zone_indices`. Entries of arrays are
numbered from 1.

## Building

The tool is built and installed with the other contrib tools. It has its own
configure script, which Neko's configure runs with the same options, and it is
compiled with the plain Fortran compiler `FC`. It needs json-fortran built
with the same compiler; if `FC` cannot use it, configure warns and the tool is
skipped.

The tool can also be compiled on its own, from its single source file, for
instance when Neko is cross-compiled for another architecture than the
machine the tool is to run on, with a json-fortran built for that machine:

```text
gfortran -O2 checkneko.F90 `pkg-config --cflags --libs json-fortran` \
    -o checkneko
```

Built this way, the banner printed at start shows the version as `unknown`.

## Schemas

The schemas are read when the tool starts, from `doc/schemas` of the Neko
source, which `make install` copies to `share/neko/schemas` under the
installation prefix. The tool finds them there, next to the `bin` directory
it is installed in, also when it is run through a symbolic link, as in a
Spack view. When it is run from the build directory of a Neko source tree, it
uses `doc/schemas` of that tree. Otherwise, for instance for a copy compiled
on its own, the directory is given with `-s`. The directory in use is printed
at start.

The validator implements the JSON Schema keywords the schemas use. When a
schema starts to use a keyword it does not implement, `checkneko` prints a
warning naming it whenever it starts, and the integration test fails, rather
than the keyword being ignored silently.
