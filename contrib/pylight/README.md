```text
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/
```

# Neko light tools (Python) {#pylight}

Standalone, CPU-only mesh utilities for Neko, in scientific Python.  The
shared library (`nekolight/`) holds the byte layout, topology table and
algorithm exactly once; five thin command-line tools sit on top of it.  No
Neko build, no MPI, no Fortran compiler; Python3 with numpy is the only
hard requirement, with scipy for partitioning (matplotlib, pyvista, ... optional).

Currently, only the following tools have python implmentations. 
A few more functionalities will be added in the near future. 

| Tool | What it does | Needs |
|---|---|---|
| `rea2nbin.py` | NEKTON `.re2` → `.nmsh` converter, byte-exact with `contrib/rea2nbin` (periodic BCs, curved elements) | numpy |
| `genmeshbox.py` | Box-mesh generator, byte-exact with `contrib/genmeshbox`; O(chunk + nelx+nely+nelz) memory | numpy |
| `mesh_checker.py` | `.nmsh` validation and diagnostics: sizes with the periodic merge (`glb_mfcs`/`glb_meds`), zones, unlabelled external faces, `--jacobian`, `--write-zone-indices` | numpy |
| `prepart.py` | Mesh partitioner (spectral / METIS / geometric) writing a reordered `.nmsh` whose linear read reproduces the partition exactly | numpy + scipy (optional: pymetis, pyamg) |
| `meshview.py` | Interactive mesh viewer: external-surface extraction + PyVista, `.vtu` export for ParaView, matplotlib fallback | numpy (optional: pyvista, matplotlib) |

Run the validation suite from this directory inside a Neko checkout:
`python3 run_tests.py` (optional dependencies that are missing SKIP their
tests rather than failing them).

## Design

* **Functions over numpy arrays.**  The one shared data structure is the
  `Mesh` named tuple (the raw `.nmsh` record arrays); there are no class
  hierarchies.  A master's student with basic numpy should be able to read
  any module top to bottom.
* **All format knowledge lives in `nekolight/formats.py`** – the `.nmsh`,
  `.re2` and `.fld` byte layouts and the vertex-ordering tables
  (`FACE_RE2`, `EDGE_RE2`, `FACET_MAP`, `SIJK`), each defined exactly once
  and documented where it is defined.
* **Atomic outputs.**  Every writer refuses an output path that names one
  of its inputs, writes to a temporary file, and renames it into place only
  after the entire run (including validation) has succeeded.  A failed run
  can never truncate an input or leave a partial output behind.
* **Validate and refuse.**  Out-of-range zone/curve references, bad labels,
  non-permutation element ids and truncated sections are hard errors in
  every tool.  Nothing silently drops or rewrites a malformed record.
* **Exact Neko semantics where it matters.**  The converter reproduces
  Neko's first-appearance, bit-exact point numbering and its fixed 3-sweep
  periodic id merge (byte-exact output); the partitioners produce block
  sizes exactly matching Neko's `linear_dist_t` (`L = nelv/P`, the first
  `nelv mod P` ranks get `L+1`), so rank *i* of a *P*-rank run receives
  exactly partition *i* – including when `nelv mod P /= 0`.
* **int64 where 8·nelv can overflow, uint64 for packed sort keys.**

## The partitioner {#pylight-prepart}

`prepart.py mesh.nmsh P [out.nmsh] [--backend spectral|metis|geometric]`

All backends share the periodic-merged, shared-vertex-weighted element dual
graph (built as one sparse E·Eᵀ product) and write contiguous per-partition
blocks in Neko's exact linear-read sizes:

* `spectral` (default): recursive spectral bisection – equivalent to Nek5000 `genmap`'s
  algorithm but with scipy eigensolvers (dense / shift-invert Lanczos /
  pyamg-preconditioned LOBPCG by sub-problem size).  Deterministic.
* `metis`: multilevel k-way via pymetis, then an exact-size repair pass
  that moves a few boundary elements to match the linear shares (the moves
  are reported). 
* `geometric`: recursive coordinate bisection on element centroids.  With
  `--no-stats` it needs neither scipy nor the vertex merge – the numpy-only
  fast path for very large meshes.  *Ignores periodic wrap-around.*

If a balanced connected bisection does not exist (e.g. a star-shaped dual
graph), one side of that bisection stays disconnected and a note is printed;
the output remains a valid permutation – only cut quality is affected.

## The viewer {#pylight-meshview}

`meshview.py mesh.nmsh [--color zone|partition|jacobian] [--nparts P]
[--export skin.vtu] [--screenshot out.png] [--matplotlib]`

The volume is never rendered: the external surface (skin) is extracted in
numpy – an n-element mesh has O(n^(2/3)) boundary quads, so even a
100M-element mesh reduces to about a million quads.  PyVista (VTK,
ParaView's rendering engine; `pip install pyvista`) then gives smooth
camera interaction, cell picking and clipping widgets; since VTK ≥ 9.4 the
same wheel renders headless (`--screenshot` on cluster nodes).  Without
pyvista, `--export skin.vtu` writes a file ParaView opens directly, and
`--matplotlib` handles small meshes (≲ 10⁴ elements) with no further
dependencies.
*Here, the vtk output is the only user tested path.*

## Memory

The whole-mesh readers hold the raw records in memory: 228 B per element
for a `.nmsh` (so ~2.3 GB at 10⁷ elements; a fat node at 10⁸+).  The two
sort-based passes are the peak consumers: point de-duplication costs about
45 B per corner (≈ 360 B/element transient) and face matching about
160 B/element.  `iter_nmsh_elements` streams the element section in chunks
for reductions that do not need the whole mesh.  `genmeshbox.py` streams
its output and needs only the three 1-D grid-line arrays plus one chunk
(a 75M×1×1 box still stores 75M+1 grid lines ≈ 600 MB – degenerate boxes
are the worst case).  Curved meshes: a `.nmsh` curve record is 532 B per
curved element; a fully curved 10⁸-element mesh carries ~53 GB of curve
records in the file itself.

## Scope

3D hex meshes only.  The Jacobian check is exact for straight-sided
elements (curve records are not applied to the geometry).  Legacy zone
types (1-4, pre-labelled-zone Neko) are carried through verbatim by the
partitioner and treated as documented boundaries by the checker; Neko's
current reader ignores them.  All tools are validated against the meshes
shipped with Neko (`run_tests.py`), but it remains your responsibility to
confirm the result is correct for your own mesh.
