# Input-output {#io}

\tableofcontents

## Mesh

Neko has it's own mesh format, `.nmsh`. All meshes should be 3D and consist of
hexahedral elements. A 2D or 1D case can be mimicked by having a single element
across selected axes and applying periodic boundary conditions.

A native mesh generator for simple box meshes, called `genmeshbox`, is part of
Neko, and is used in some of the examples, for example `advecting_cone`. The
usage is straightforward, and is well described by the help string provided when
running the utility.

The main way of obtaining a `.nmsh` is converting a Nek5000 `.re2` mesh file
using the `rea2nbin` utility. Nek5000, in turn, has several converters among its
tools, such as `gmsh2nek`. The latter is also available under the `contrib`
directory in Neko. The workflow is thus to export your mesh into a format, which
can be converted to `.re2`, and then convert the `.re2` to `.nmsh`. Of course
one can also use native Nek5000 tools to produce the `.re2`, such as `genbox`.
In the future, native support for other formats than `.nmsh` will be added for
convenience.

It should be noted, that the `.re2` format allows to store boundary conditions.
This is relevant for users that have old Nek5000 cases, who wish to convert them
to Neko. Boundary condition is converted and used by Neko, and for these
boundaries one does not need to provide information in the `boundary_types`
keyword in the case file. This is why in some of the examples, `boundary_types`
is only filled for some of the boundaries. However, since this feature
complicates the code and leads to somewhat confusing case setups, it is planned
to deprecate it at some point. Note that periodic boundaries are also directly
encoded into the mesh file, and this will remain so in the future.

## Three-dimensional field output
Neko stores the 3D fields with results in the `.fld` format, which is the same
as in Nek5000. The advantage of adopting this format, is that there is a reader
in Paraview and Visit, which can be used to visualize them. Note that the latest
version of Paraview actually has two reader for `.fld`. For now, Neko has only
been tested with the older reader, which uses Visit under the hood.  A file with
the `.nek5000` extension is used as the entry point for the readers and stores
some metadata. Users may also find the Python package `pymech` useful for
working with `.fld`s. Note that only the first output `.fld` file stores the
mesh.

## Checkpoint files
Simulations cannot be restarted from `.fld` files. Instead, separate checkpoint
files can be output for the purpose of restarts. These contain additional
information allowing a clean restart, with, e.g., the correct time integration
order. A separate file format, `.chkp` is adopted for the checkpoint files.
