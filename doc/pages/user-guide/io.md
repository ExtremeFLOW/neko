# Input-output {#io}

\tableofcontents

## Mesh

Neko has its own mesh format, `.nmsh`. All meshes should be 3D and consist of
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
to Neko. Boundary conditions are not converted by Neko as such, but the faces
where they were applied are saved into regular boundary zones. The relevant
conditions should then be prescribed in the case file as usual. Note also that
periodic boundaries are directly encoded into the mesh file, and this will
remain so in the future.

## Three-dimensional field output
Neko stores the 3D fields with results in the `.f#####` format, which is the same
as in Nek5000. The advantage of adopting this format, is that there is a reader
in Paraview and VisIt, which can be used to visualize them. Note that the latest
version of Paraview actually has two readers for `.f#####`. For now, Neko has
been thoroughly tested with the older reader, which uses VisIt under the hood.
However, the new reader was used for the smaller cases as well. A file with
the `.nek5000` extension is used as the entry point for the readers and stores
some metadata. Users may also find the Python package
[`pysemtools`](https://github.com/ExtremeFLOW/pySEMTools) useful for
working with `.f#####`s. Note that only the first output `.f#####` file stores the
mesh. Detailded description of the file format can be found in
[Nek5000 documentation](https://nek5000.github.io/NekDoc/problem_setup/case_files.html#restart-output-files-f).
Notice, for brievity in this manual we often call `.f#####` file format an `.fld` one
and in some cases Neko even outputs files with `.fld` extension. However, in all cases
we refer to the binary `.f#####` format and not to the Nek5000 text file format with
the same extension.

### Compression of field output
Neko supports compression of the 3D field data when writing through the ADIOS2.
If the ADIOS2 dependency was compiled with compression library support like
BigWhoop, SZ, or ZFP the field output can be compressed using one of these lossy
compressors (and other lossless compression techniques.) Using algorithms based
on image compression, such as BigWhoop and ZFP, the local structure from higher
order elements needs to be regarded. Setting `output\_layout=2` uses the
implicit neighbourhood relations of DoFs in x,y,z in a 4D layout. The first
three dimensions represent the DoFs and the fourth dimension refers to the
number of elements in numeric order. This layout can be used for BigWhoop and
ZFP. Additionally, benefits in compression efficiency can be observed for
BigWhoop when packing and compressing field parameters together, effectively
adding a fifth dimension for the number of field parameters. This layout is
enabled when `output\_format=3`. (The default `output_layout=1` creates a
contiguous 1D array access pattern equivalent to the data layout in the nek5000
files.)

The compression algorithms are controlled by an additional file `adios2.xml`.
Please refer to the documentation of ADIOS2 (and the specific compression
library) to configure the data "operators" of the parallel I/O library.

The tool `adios2_to_nek5000` can be used to decompress the field data for
postprocessing and visualization with the conventional nek5000 format.
Additionally, the data can be compressed using this tool as a postprocessing
step using the `adios2.xml` configuration and a given uncompressed data set as
input.
```
adios2_to_nek5000 input.bp/fld output.fld/bp .false. 1
```
The bool parameter specifies whether output is written using double precision.
The parsed integer specifies the data layout in the case of an ADIOS2 `.bp`
output file.

The tool `psnr` can be used to analyze the Peak-Signal-to-Noise Ratio (PSNR),
quantifying the compression efficiency, namely, the relation between compression
ratio and the lost accuracy due to lossy data compression.
```
psnr compressed_fields.bp uncompressed_fields.fld/bp .false.
```
(The bool parameter specifies whether output should be double precision.) Test
data sets can be generated during preprocessing and in an initial case setup
workflow to compare and find an optimal compression rate with an acceptable
accuracy loss. The value of the PSNR decreases with increasing accuracy loss.
Tune the compression by increasing the compression ratio until a desired lower
limiting value for the PSNR is reached. As a rule of thumb, target values of 60
or 40 can be used as lower limit for accurate postprocessing or visualization,
respectively.  Separately, it is recommended to check for compression errors
from lossy compressors in the specific quantities of interest in postprocessing
and visualization.

## Checkpoint files
Simulations cannot be restarted from `.fld` files (although you can use an `fld`
to provide initial conditions). Instead, separate checkpoint files can be output
for the purpose of restarts. These contain additional information allowing a
clean restart, with, e.g., the correct time integration order. A separate file
format, `.chkp` is adopted for the checkpoint files.
