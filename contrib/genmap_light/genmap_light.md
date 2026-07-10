```text
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/
```

# genmap_light

A standalone, CPU-only **spectral-bisection mesh partitioner** for Neko `.nmsh`
files — a lightweight reimplementation of Nek5000's `genmap` that reads an
`.nmsh`, partitions it into `nparts` parts, and writes a reordered `.nmsh` whose
elements are grouped into contiguous per-partition blocks, exactly like
`contrib/prepart` (so a later *linear* read reproduces the partition). Curves and
zones are first-class.

Needs only a Fortran compiler (no MPI, no Neko library, no ParMETIS):

```sh
make                                          # gfortran -O2 -o genmap_light genmap_light.f90
./genmap_light mesh.nmsh 8                     # -> mesh_8.nmsh (8 parts)
./genmap_light mesh.nmsh 16 out.nmsh
make check                                     # partition + reorder validation vs a SciPy oracle
```

Usage mirrors `prepart`:

```
genmap_light mesh.nmsh nparts [out.nmsh]
```

`nparts` is any positive integer (need not be a power of two); the default output
is `<base>_<nparts>.nmsh`.

> `make check` reads meshes from a Neko checkout and needs `numpy`+`scipy`, so it
> expects this folder inside a Neko source tree (`../../examples`, `../../tests`);
> otherwise pass the tree explicitly:
> `python3 validate_genmap.py ./genmap_light /path/to/neko [nparts]`.
> Partitioning your own mesh needs no Neko tree and no Python.

## Why it exists

Neko partitions meshes with **ParMETIS** (via `prepart`, and the runtime load
balancer). ParMETIS is excellent, but it is an external MPI library you have to
build and link. Nek5000's `genmap` uses a lighter-weight, self-contained
algorithm — **recursive spectral bisection (RSB)** using the Fiedler vector of
the element graph — that is often nearly as good and needs no dependencies.
`genmap_light` brings that algorithm to Neko's `.nmsh` format as a small
standalone tool.

Like the other `*_light` tools it does **no Neko initialisation** and links
nothing, so it runs anywhere with just a Fortran compiler — in particular you do
not need a separate CPU build of a GPU-configured Neko, and you do not need
ParMETIS/MPI available just to (re)partition a mesh offline.

> **Tested, but please check your own case.** The partition quality is
> cross-checked against an independent SciPy exact-Fiedler eigensolver and the
> reorder is checked to be a pure permutation (geometry, curves and zones
> preserved) across the meshes shipped with Neko (see [Validation](#validation)).
> It remains your responsibility to confirm the result is correct for *your* mesh
> — `make check` on a case you trust first.

## The algorithm (faithful to genmap)

1. **Element dual graph.** Each element is a graph node; two elements are
   connected if they share one or more vertices, with **edge weight = number of
   shared vertices** (face-neighbours share 4, edge-neighbours 2, corner-
   neighbours 1) — exactly Nek5000 genmap's `c2c` weighting. This becomes the
   weighted graph Laplacian `L` (`L_ii = Σ weights`, `L_ij = −weight`).
   Point ids come straight from the `.nmsh`; **periodic** facets are merged first
   using the stored `glb_pt_ids` (the same merge Neko applies on read), so
   periodic neighbours are graph-adjacent.
2. **Low spectral modes.** The few lowest non-trivial eigenvectors of `L`
   (the Fiedler vector and the next few) — computed by a Lanczos iteration on `L`
   that deflates the constant null vector, just as genmap does. `genmap_light`
   uses **full reorthogonalisation** and a generic (all-mode) start, which is more
   robust than genmap's plain Lanczos: plain Lanczos loses orthogonality and can
   under-resolve the low modes on high-aspect-ratio or periodic meshes and then
   bisect along the wrong one.
3. **Bisection (best of the low modes).** For each of those lowest modes, sort
   the elements by the mode value, split at the balance point, and measure the
   resulting edge cut; keep the **lowest-cut** split, preferring splits whose two
   halves are both connected (a DFS check). If none is connected, a
   region-growing (BFS) bisection restores connectivity (genmap's
   connectivity-preserving strategy). Trying the few lowest modes — not just the
   single Fiedler vector — makes the bisection robust on meshes with a
   degenerate/near-degenerate Fiedler space (cylinders, symmetric boxes), where
   the specific eigenvector chosen determines the cut.
4. **Recurse** on each half until `nparts` parts remain. Non-power-of-two
   `nparts` is handled by splitting the target part-count proportionally, so all
   final parts differ by at most one element.
5. **Reorder + write** (the `prepart` contract): number the parts' elements into
   contiguous blocks, write a new `.nmsh` with `el_idx` renumbered `1..nelv`,
   point ids and coordinates **unchanged**, and every curve and zone carried
   through with its element reference remapped. A later linear read on `nparts`
   ranks then lands each part on its rank.

**Curves are first-class.** Every curve record is carried through: the element id
is remapped, and the geometric payload (`curve_data(5,12)`, `curve_type(12)`) is
copied **verbatim** — coordinates are never transformed, curvature never dropped.
Periodic (type-5) zones keep their `glb_pt_ids` verbatim (those reference point
ids, which the reorder does not change) with only the element references
remapped; labeled (type-7) zones keep their label.

## Validation

`validate_genmap.py` runs genmap_light on every 3D mesh in a Neko tree and checks
two things.

**1. The reorder is correct (the prepart contract).** For every mesh the output
is a **pure permutation** of the input — element geometry, all curve payloads,
and all zone records are preserved as multisets — with `el_idx` renumbered to a
contiguous `1..nelv` and the elements grouped into contiguous per-partition
blocks, so a linear read reproduces the partition. Curves and periodic/labeled
zones are verified carried through with only element references remapped.

**2. The partition is good.** It rebuilds the same weighted dual-graph Laplacian
and computes an **independent recursive spectral bisection with SciPy's exact
eigensolver** (`scipy.sparse.linalg.eigsh`), then compares genmap_light's
weighted edge cut and balance against that oracle and against a naive linear
split. genmap_light's cut matches the exact-Fiedler oracle closely (often equal
or better, thanks to the connectivity-preserving fallback) and is far below the
linear baseline where spectral partitioning helps, with balance exact
(part sizes differ by ≤ 1).

Representative results (`nparts = 8`; weighted shared-vertex edge cut):

| mesh | nelv | genmap_light | SciPy exact-Fiedler | linear |
|---|--:|--:|--:|--:|
| `turb_pipe` (curved+periodic) | 36480 | **49100** | 50336 | 199120 |
| `tgv/32768` (3× periodic) | 32768 | **111181** | 101954 | 131072 |
| `turb_channel` | 5832 | **25886** | 26624 | 41216 |
| `cyl_boundary_layer` | 6939 | **13376** | 13448 | 29262 |
| `cylinder` | 1824 | **6192** | 6192 | 8088 |
| `hemi` | 2042 | **4900** | 4960 | 24380 |
| `poisson/16384` | 16384 | **36456** | 31868 | 111132 |

Across all 42 shipped 3D meshes the reorder is a pure permutation of the input,
balance is exact (part sizes differ by ≤ 1), and the cut is at or near the exact
spectral optimum (often below it, thanks to the best-of-modes choice). Against a
naive linear split the cut is typically 1.2–5× lower, up to ~9× on stretched
meshes, and equal only where the file order is already optimal (extruded,
1-D-like boxes).

## Scope / limitations

- **3D hex meshes** (genmap's common case). 2D (`gdim=2`) exits with a clear
  error.
- **Serial / offline tool** (like `prepart`): it loads the whole mesh, so peak
  memory is roughly 2 kB per element (element records + the weighted dual graph
  + the Lanczos vectors) — small next to a solver run, but it is not itself
  distributed. Practical up to meshes of tens of millions of elements; the CSR
  dual-graph index is 32-bit, as in genmap.
- **Partition quality vs ParMETIS.** RSB gives a good, low-cut partition but
  ParMETIS's multilevel k-way method can still edge it out on some meshes; use
  `prepart` when you have ParMETIS and want the last few percent. genmap_light is
  the dependency-free option.
- **Python sibling.** `prepart_light` (shipped alongside) implements the same
  reorder contract in Python with three backends (spectral / METIS / geometric);
  it is usually faster and, with `pymetis`, gives the best cut. Use
  `genmap_light` when a Fortran compiler is all you have, `prepart_light` when
  Python with NumPy/SciPy is available.
- **Determinism.** The Fiedler start is a fixed deterministic hash, so a given
  mesh + `nparts` always yields the same output.
- **Diagnostics.** Setting the `GENMAP_DEBUG` environment variable prints one
  line per bisection to stderr (subset size, split point, connectivity, cut).
- **Legacy zone types.** Zone records of types other than 5/7 (found in some
  older meshes) are carried through verbatim with the element id renumbered.
- The only non-geometric fields that differ from a `prepart` run are ordering and
  the labeled-zone `p_e`/`glb_pt_ids` that Neko leaves uninitialised (written as
  deterministic zeros here); the reader ignores them.

## Files

- `genmap_light.f90` — the partitioner (hash + RSB numerics modules + program).
- `validate_genmap.py` — the SciPy-oracle + permutation/curve/zone harness
  (`make check`; needs `numpy`+`scipy`).
- `Makefile` — `make`, `make check`, `make clean`.
