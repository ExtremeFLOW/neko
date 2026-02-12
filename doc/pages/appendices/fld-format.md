# Neko .fld field format {#fld-format}

\tableofcontents

This appendix documents the Nek5000 file format as written and read by Neko in
``src/io/fld_file.f90``. We refer to this format as `.fld` for brevity. Neko
writes ``.fld`` files using MPI-IO, and the layout on disk follows the code
paths in that module.

## File naming and series

Neko writes a file series in Nek5000 style:

- Data files: ``<base>0.f#####`` (zero-padded index)
- Metadata file: ``<base><start_counter>.nek5000``

The ``.nek5000`` file contains the file template and time range for the series
and is used by Neko when reading a series.

## File overview

Each ``.fld`` data file is binary and consists of:

1. **Header** (132 characters)
2. **Test pattern** (single precision, value ``6.54321``)
3. **Element index list** (``glb_nelv`` integers)
4. **Field blocks** (mesh, velocity, pressure, temperature, scalars)
5. **Metadata blocks** (per-element min/max, only for 3D)

The presence and ordering of the field blocks is described by ``rdcode`` in the
header.

## Header

The header is a 132-character record written with the format:

```
#std <fld_data_size> <lx> <ly> <lz> <glb_nelv> <glb_nelv> <time> <step>
     <1> <1> <rdcode[10]>
```

In code, this is:

```
'#std', i1, i2, i2, i2, i10, i10, e20.13, i9, i6, i6, 10a
```

Fields:

- ``fld_data_size``: byte size of field data (``MPI_REAL_SIZE`` or
  ``MPI_DOUBLE_PRECISION_SIZE``). In practice, this is either 4 or 8.
- ``lx, ly, lz``: polynomial orders (GLL grid size per element direction)
- ``glb_nelv``: global number of elements (written twice)
- ``time``: simulation time
- ``step``: output counter
- The two ``i6`` values are currently written as ``1``
- ``rdcode``: 10-character field descriptor

### rdcode

``rdcode`` is a compact descriptor of which fields are present and in which
order:

- ``X``: mesh coordinates
- ``U``: velocity (vector)
- ``P``: pressure (scalar)
- ``T``: temperature (scalar)
- ``Snn``: scalars, with a two-digit count (``S00`` .. ``S99``)

The writer appends ``S`` followed by two digits for the scalar count. The reader
parses these digits to determine the number of scalar fields.

## Test pattern

Immediately after the header, the writer stores a single-precision real with
value ``6.54321``. The reader validates this value to guard against format or
endianness mismatches.

## Element index list

After the test pattern, the file contains ``glb_nelv`` integers storing the
global element ids. Each MPI rank reads/writes its own slice based on
``offset_el``.

## Field blocks

Field data is written in the order indicated by ``rdcode``. Within each block,
Neko writes data in **element-major** order:

- Loop elements ``e = 1..nelv``
- For each element, loop GLL points ``j = 1..lxyz``

### Vector fields (X and U)

For vector fields, the layout per element is:

1. All ``x`` values (``j = 1..lxyz``)
2. All ``y`` values (``j = 1..lxyz``)
3. All ``z`` values (``j = 1..lxyz``) if ``gdim == 3``

The block size is ``gdim * lxyz * glb_nelv`` values. Note that if `U` is
written, all three components have to be present. One cannot, for example, write
only one component of velocity.

If `X` is not written, then the files does not contain the mesh!
### Scalar fields (P, T, S)

Scalar fields are stored as ``lxyz * glb_nelv`` values in element-major order.
Scalars are written in sequence ``S1, S2, ...``.

### Precision

Field blocks are written in either single or double precision, depending on the
``fld_data_size`` header value. Neko writes single precision by default unless
configured otherwise.

## Metadata blocks (3D only)

For 3D output (``gdim == 3``), Neko appends per-element min/max metadata for
all written fields. These metadata blocks are always **single precision** and
follow the same field ordering as the data blocks.

- Vector fields: ``2 * gdim * glb_nelv`` values
  - ``xmin, xmax, ymin, ymax, zmin, zmax`` per element
- Scalar fields: ``2 * glb_nelv`` values
  - ``min, max`` per element

For 2D output, these metadata blocks are not written.

## 2D fields and Z-velocity

When writing 2D fields, Neko sets ``gdim = 2`` (``lz = 1``). If a Z-velocity
component exists, it is stored as an additional scalar field in the output.
This is a semantic convention rather than a special file layout.

## Portability notes

The ``.fld`` format is not self-describing beyond the header. It depends on:

- Endianness (typically little).
- ``MPI_INTEGER`` size
- ``MPI_REAL`` and ``MPI_DOUBLE_PRECISION`` sizes

For robust I/O, use Neko's ``fld_file`` reader/writer or convert using tools
that understand the Nek5000 format.
