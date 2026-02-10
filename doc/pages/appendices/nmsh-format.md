# Neko .nmsh mesh format {#nmsh-format}

\tableofcontents

This appendix documents the Neko binary mesh format (``.nmsh``) as implemented
in ``src/io/format/nmsh.f90`` and ``src/io/nmsh_file.f90``. The format is
written and read with MPI-IO using MPI derived types, so the on-disk layout
mirrors the in-memory layout of the Fortran derived types.

## File overview

An ``.nmsh`` file is a binary stream with the following blocks, in order:

1. **Header** (two ``MPI_INTEGER`` values)
   - ``nelv``: number of elements (global)
   - ``gdim``: geometric dimension (2 or 3)

2. **Element block**
   - ``nelv`` records of ``nmsh_quad_t`` if ``gdim == 2``
   - ``nelv`` records of ``nmsh_hex_t`` if ``gdim == 3``

3. **Zone block**
   - ``nzones``: number of zone records (``MPI_INTEGER``)
   - ``nzones`` records of ``nmsh_zone_t``

4. **Curve block**
   - ``ncurves``: number of curve records (``MPI_INTEGER``)
   - ``ncurves`` records of ``nmsh_curve_el_t``

If ``nzones`` or ``ncurves`` is zero, the corresponding record array is omitted.

## Data types and portability

The file uses:

- ``MPI_INTEGER`` for integer fields
- ``MPI_DOUBLE_PRECISION`` for real fields (``real(kind=dp)``)

The layout is **not self-describing**. It depends on the platform endianness,
``MPI_INTEGER`` size, and Fortran derived-type padding/alignment. In practice,
all ``.nmsh`` files in this repository are little-endian with 32-bit integers,
but external readers should **not** assume this. For robust I/O, use Neko's
``nmsh_file`` reader/writer or convert via ``rea2nbin``.

All element and vertex indices stored in the file are **global, 1-based** ids
as used by Neko internally.

## Element records

### Vertex record (``nmsh_vertex_t``)

| Field | Type | Description |
| --- | --- | --- |
| ``v_idx`` | ``MPI_INTEGER`` | Global vertex id |
| ``v_xyz`` | ``MPI_DOUBLE_PRECISION[3]`` | Vertex coordinates (x, y, z) |

### Quad record (``nmsh_quad_t``)

| Field | Type | Description |
| --- | --- | --- |
| ``el_idx`` | ``MPI_INTEGER`` | Global element id |
| ``v(1:4)`` | 4×``nmsh_vertex_t`` | Quad vertices |

### Hex record (``nmsh_hex_t``)

| Field | Type | Description |
| --- | --- | --- |
| ``el_idx`` | ``MPI_INTEGER`` | Global element id |
| ``v(1:8)`` | 8×``nmsh_vertex_t`` | Hex vertices |

### Vertex ordering

The writer stores vertices using the ``vcyc_to_sym`` mapping in
``src/io/nmsh_file.f90``:

```
[1, 2, 4, 3, 5, 6, 8, 7]
```

This is the ordering used for the ``nmsh_hex_t`` / ``nmsh_quad_t`` records. The
reader swaps vertices when constructing elements to restore Neko's internal
ordering. External tools should preserve the file ordering as-is.

## Zone records (boundary and periodic data)

The zone record is defined by ``nmsh_zone_t``:

| Field | Type | Description |
| --- | --- | --- |
| ``e`` | ``MPI_INTEGER`` | Global element id |
| ``f`` | ``MPI_INTEGER`` | Facet number (1-based) |
| ``p_e`` | ``MPI_INTEGER`` | Periodic partner element id |
| ``p_f`` | ``MPI_INTEGER`` | Periodic partner facet id (or label id) |
| ``glb_pt_ids(1:4)`` | ``MPI_INTEGER[4]`` | Global point ids for the facet |
| ``type`` | ``MPI_INTEGER`` | Zone type |

Currently, Neko writes:

- ``type = 5``: periodic facet. ``p_e`` and ``p_f`` identify the partner facet;
  ``glb_pt_ids`` stores the ordered global point ids of the facet.
- ``type = 7``: labeled facet. ``p_f`` stores the zone label index.

Other zone types are not emitted by the current writer.

## Curve records (curved edges/faces)

The curve record is defined by ``nmsh_curve_el_t``:

| Field | Type | Description |
| --- | --- | --- |
| ``e`` | ``MPI_INTEGER`` | Global element id |
| ``curve_data`` | ``MPI_DOUBLE_PRECISION[5,12]`` | Curve parameters |
| ``type`` | ``MPI_INTEGER[12]`` | Curve type per edge |

The meaning of ``curve_data`` and the ``type`` codes follows the same
conventions as the ``re2`` importer in ``src/io/re2_file.f90``. The current
implementation recognizes:

- ``0``: no curve
- ``3``: arc (``'C'`` in ``.re2``)
- ``4``: mapped curve (``'m'`` in ``.re2``)

Other curve types (e.g. ``1`` / ``2`` from ``'s'`` / ``'e'``) are read but
treated as unsupported by the importer.

## 2D meshes

If the header ``gdim`` is ``2``, the file stores ``nmsh_quad_t`` records. The
reader in ``nmsh_file_read_2d`` extrudes these into a thin 3D slab:

- Bottom vertices are placed at ``z = 0``
- Top vertices are placed at ``z = 1``
- Periodic facets are added between the two planes

This behavior is specific to the reader and is **not** encoded in the file.
